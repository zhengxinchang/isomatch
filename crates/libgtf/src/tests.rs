use std::{
    fs::{self, File},
    io::Cursor,
    time::{SystemTime, UNIX_EPOCH},
};

use crate::gtf::Strand;
use crate::index::tx::TxBaseLoadArgs;
use crate::index::{
    AttrIndexReader, ChromDirectoryEntry, ISOMX_VERSION, IndexHeader, IndexReader, JunctionPool,
    JunctionSpan, SpliceSitePool, SpliceSiteSpan, StringPool, StringSpan, TxBase,
};
use crate::traits::{Decodable, DiskSize, Encodable, PartialLoad};
use crate::{BuildConfig, BuildEvent, build_index_with_events};

#[test]
fn tx_base_disk_layout_round_trip() {
    let tx = TxBase::new(
        0x0102_0304_0506_0708,
        7,
        100,
        350,
        Strand::Minus,
        0x0102_0304_0506_0708_090a_0b0c_0d0e_0f10,
        0x1020_3040_5060_7080_90a0_b0c0_d0e0_f000,
        3,
        SpliceSiteSpan {
            offset: 11,
            count: 2,
        },
        JunctionSpan {
            offset: 23,
            count: 4,
        },
        StringSpan {
            offset: 37,
            byte_len: 5,
        },
        StringSpan {
            offset: 53,
            byte_len: 6,
        },
    )
    .unwrap();

    let bytes = tx.encode().unwrap();
    assert_eq!(TxBase::DISK_SIZE, 96);
    assert_eq!(bytes.len(), TxBase::DISK_SIZE);
    assert_eq!(&bytes[0..8], &0x0102_0304_0506_0708u64.to_le_bytes());
    assert_eq!(&bytes[8..12], &100u32.to_le_bytes());
    assert_eq!(&bytes[12..16], &350u32.to_le_bytes());
    assert_eq!(&bytes[16..18], &5u16.to_le_bytes());
    assert_eq!(&bytes[52..56], &23u32.to_le_bytes());
    assert_eq!(&bytes[56..58], &4u16.to_le_bytes());
    assert_eq!(&bytes[58..62], &11u32.to_le_bytes());
    assert_eq!(&bytes[62..64], &2u16.to_le_bytes());
    assert_eq!(&bytes[64..72], &37u64.to_le_bytes());
    assert_eq!(&bytes[72..80], &5u64.to_le_bytes());
    assert_eq!(&bytes[80..88], &53u64.to_le_bytes());
    assert_eq!(&bytes[88..96], &6u64.to_le_bytes());

    let loaded = TxBase::load_range(
        &mut Cursor::new(bytes.as_slice()),
        0,
        TxBase::DISK_SIZE,
        TxBaseLoadArgs { chrom_id: 7 },
    )
    .unwrap();
    assert_eq!(loaded.encode().unwrap(), bytes);
    assert_eq!(loaded.tx_idx(), tx.tx_idx());
    assert_eq!(loaded.chrom_id(), tx.chrom_id());
    assert_eq!(loaded.tx_boundary(), tx.tx_boundary());
}

#[test]
fn index_header_disk_layout_round_trip() {
    let mut header = IndexHeader::new(2, 1234, 4096, [0xab; 16], true, false, 12, 3, 8);
    header.total_tx_n = 42;
    header.reserved_to_4k[0] = 0x5a;
    header.reserved_to_4k[header.reserved_to_4k.len() - 1] = 0xa5;

    let bytes = header.encode().unwrap();
    assert_eq!(IndexHeader::DISK_SIZE, 4096);
    assert_eq!(bytes.len(), IndexHeader::DISK_SIZE);
    assert_eq!(&bytes[0..4], b"ISOM");
    assert_eq!(&bytes[4..8], &ISOMX_VERSION.to_le_bytes());
    assert_eq!(&bytes[56..64], &header.total_tx_n.to_le_bytes());
    assert_eq!(bytes[72], 0x5a);
    assert_eq!(bytes[4095], 0xa5);

    let loaded = IndexHeader::decode_from(&mut Cursor::new(bytes), ()).unwrap();
    assert_eq!(loaded.magic, header.magic);
    assert_eq!(loaded.version, header.version);
    assert_eq!(loaded.flags.bits, header.flags.bits);
    assert_eq!(loaded.chrom_count, header.chrom_count);
    assert_eq!(loaded.gtf_file_size, header.gtf_file_size);
    assert_eq!(loaded.index_file_size, header.index_file_size);
    assert_eq!(loaded.md5, header.md5);
    assert_eq!(loaded.chrom_name_table_len, header.chrom_name_table_len);
    assert_eq!(loaded.total_tx_n, header.total_tx_n);
    assert_eq!(loaded.missing_seqid_count, header.missing_seqid_count);
    assert_eq!(
        loaded.missing_seqid_table_len,
        header.missing_seqid_table_len
    );
    assert_eq!(loaded.reserved_to_4k, header.reserved_to_4k);
}

#[test]
fn chrom_directory_entry_disk_layout_round_trip() {
    let entry = ChromDirectoryEntry {
        chrom_id: 7,
        chrom_name_offset: 11,
        chrom_name_len: 13,
        global_tx_offset: 17,
        global_tx_count: 19,
        global_junction_pool_offset: 23,
        global_junction_pool_len: 29,
        global_string_pool_offset: 31,
        global_string_pool_len: 37,
        global_splice_site_pool_offset: 41,
        global_splice_site_pool_len: 43,
    };

    let bytes = entry.encode().unwrap();
    assert_eq!(ChromDirectoryEntry::DISK_SIZE, 74);
    assert_eq!(bytes.len(), ChromDirectoryEntry::DISK_SIZE);
    assert_eq!(&bytes[0..2], &entry.chrom_id.to_le_bytes());
    assert_eq!(&bytes[2..6], &entry.chrom_name_offset.to_le_bytes());
    assert_eq!(&bytes[6..10], &entry.chrom_name_len.to_le_bytes());
    assert_eq!(&bytes[10..18], &entry.global_tx_offset.to_le_bytes());
    assert_eq!(
        &bytes[66..74],
        &entry.global_splice_site_pool_len.to_le_bytes()
    );

    let loaded = ChromDirectoryEntry::decode_from(&mut Cursor::new(bytes), ()).unwrap();
    assert_eq!(loaded.chrom_id, entry.chrom_id);
    assert_eq!(loaded.chrom_name_offset, entry.chrom_name_offset);
    assert_eq!(loaded.chrom_name_len, entry.chrom_name_len);
    assert_eq!(loaded.global_tx_offset, entry.global_tx_offset);
    assert_eq!(loaded.global_tx_count, entry.global_tx_count);
    assert_eq!(
        loaded.global_junction_pool_offset,
        entry.global_junction_pool_offset
    );
    assert_eq!(
        loaded.global_junction_pool_len,
        entry.global_junction_pool_len
    );
    assert_eq!(
        loaded.global_string_pool_offset,
        entry.global_string_pool_offset
    );
    assert_eq!(loaded.global_string_pool_len, entry.global_string_pool_len);
    assert_eq!(
        loaded.global_splice_site_pool_offset,
        entry.global_splice_site_pool_offset
    );
    assert_eq!(
        loaded.global_splice_site_pool_len,
        entry.global_splice_site_pool_len
    );
}

#[test]
fn junction_pool_encode_load_and_lookup() {
    let mut pool = JunctionPool::new();
    let first = pool.add(&[150, 200]).unwrap();
    let second = pool.add(&[250, 300, 350, 400]).unwrap();
    let bytes = pool.encode().unwrap();
    assert_eq!(bytes.len(), 6 * 4);
    assert_eq!(&bytes[0..4], &150u32.to_le_bytes());

    let loaded =
        JunctionPool::load_range(&mut Cursor::new(bytes.clone()), 0, bytes.len(), 7).unwrap();
    assert_eq!(loaded.get(first).unwrap(), &[150, 200]);
    assert_eq!(loaded.get(second).unwrap(), &[250, 300, 350, 400]);
}

#[test]
fn string_pool_encode_load_and_lookup() {
    let mut pool = StringPool::new();
    let tx_id = pool.add("tx-1").unwrap();
    let gene_id = pool.add("g\u{e8}ne-1").unwrap();
    assert_eq!(pool.add("tx-1").unwrap(), tx_id);

    let bytes = pool.encode().unwrap();
    assert_eq!(bytes, "tx-1g\u{e8}ne-1".as_bytes());
    let loaded =
        StringPool::load_range(&mut Cursor::new(bytes.clone()), 0, bytes.len(), ()).unwrap();
    assert_eq!(loaded.get(tx_id).unwrap(), "tx-1");
    assert_eq!(loaded.get(gene_id).unwrap(), "g\u{e8}ne-1");
}

#[test]
fn splice_site_pool_encode_load_and_lookup() {
    let mut pool = SpliceSitePool::new();
    let first = pool
        .add_pairs(&[(b"GT".to_vec(), b"AG".to_vec())], Strand::Plus)
        .unwrap();
    let second = pool
        .add_pairs(&[(b"GC".to_vec(), b"AG".to_vec())], Strand::Plus)
        .unwrap();

    let bytes = pool.encode().unwrap();
    assert_eq!(bytes, [0x01, 0x21]);
    let loaded =
        SpliceSitePool::load_range(&mut Cursor::new(bytes.clone()), 0, bytes.len(), ()).unwrap();
    assert_eq!(
        loaded.get_pair(first).unwrap(),
        pool.get_pair(first).unwrap()
    );
    assert_eq!(
        loaded.get_pair(second).unwrap(),
        pool.get_pair(second).unwrap()
    );
    assert!(loaded.get_pair(first).unwrap()[0].is_canonical());
}

#[test]
fn build_index_pipeline_reports_events_and_writes_readable_indexes() {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let root = std::env::temp_dir().join(format!("libgtf-build-{}-{stamp}", std::process::id()));
    let temp_parent = root.join("tmp");
    fs::create_dir_all(&temp_parent).unwrap();

    let fasta_path = root.join("reference.fa");
    let sequence = "AAAAAGTAAGCCCCCCCCCCAAAAAAAAAAAAAAAAAAAA";
    fs::write(&fasta_path, format!(">chr1\n{sequence}\n")).unwrap();
    fs::write(
        root.join("reference.fa.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        ),
    )
    .unwrap();

    let gtf_path = root.join("input.gtf");
    fs::write(
        &gtf_path,
        concat!(
            "chr1\ttest\ttranscript\t1\t20\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\ttest\texon\t1\t5\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\ttest\ttranscript\t30\t35\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
            "chr1\ttest\texon\t30\t35\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
        ),
    )
    .unwrap();

    let isomx_path = root.join("output.isomx");
    let isoms_path = root.join("output.isoms");
    let config = BuildConfig {
        gtf_path,
        reference_fasta: fasta_path,
        transcript_fasta: None,
        isomx_output: isomx_path.clone(),
        isoms_output: isoms_path.clone(),
        skip_missing_reference_seqids: false,
        temp_dir: Some(temp_parent.clone()),
    };
    let mut events = Vec::new();
    let report = build_index_with_events(&config, |event| events.push(event)).unwrap();

    assert_eq!(report.transcript_count, 2);
    assert_eq!(report.gene_count, 2);
    assert_eq!(report.plus_strand_transcript_count, 1);
    assert_eq!(report.minus_strand_transcript_count, 1);
    assert_eq!(report.mono_exon_transcript_count, 1);
    assert_eq!(report.multi_exon_transcript_count, 1);
    assert_eq!(report.all_canonical_transcript_count, 1);
    assert_eq!(report.junction_count, 1);
    assert_eq!(report.canonical_junction_count, 1);
    assert_eq!(
        events,
        vec![
            BuildEvent::LoadingFasta,
            BuildEvent::IndexingGtf,
            BuildEvent::InitializingBuilders,
            BuildEvent::ProcessingChromosome {
                name: "chr1".to_string(),
            },
            BuildEvent::Finalizing,
        ]
    );

    let index = IndexReader::open(File::open(&isomx_path).unwrap(), 0).unwrap();
    assert_eq!(index.transcript_count(), 2);
    assert_eq!(index.chromosome_names(), ["chr1"]);
    let mut attributes = AttrIndexReader::open(&isoms_path).unwrap();
    assert!(attributes.get_attr(0).unwrap().is_some());
    assert!(temp_parent.read_dir().unwrap().next().is_none());

    fs::remove_dir_all(root).unwrap();
}
