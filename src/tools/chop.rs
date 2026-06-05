use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

use flate2::{Compression, write::GzEncoder};
use log::info;
use num_format::{Locale, ToFormattedString};

use crate::{
    ChopArgs, ChopMode,
    tools::tools_error::ToolError,
    traits::ArgValidate,
    utils::{greetings2, open_file_bufread, require_file},
};

impl ArgValidate for ChopArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file("Input GTF", &self.input, &mut error_msg, &mut has_error);

        if has_error {
            panic!("{}", error_msg);
        }
    }
}

pub fn run_chop(args: &ChopArgs) -> Result<(), ToolError> {
    greetings2(&args);

    args.validate();

    info!("Loading GTF");
    let mut reader = open_file_bufread(&args.input)?;
    let mut out_path = args.out.clone();
    out_path.add_extension("chopped.gtf.gz");
    let mut writer = BufWriter::new(GzEncoder::new(
        File::create(&out_path)?,
        Compression::default(),
    ));

    let keep_attrs = parse_keep_attrs(args.keep_attrs.as_deref(), args.keep_check_case);

    info!("Start chop GTF");
    let mut line = Vec::new();
    let mut line_no = 0;
    while reader.read_until(b'\n', &mut line)? != 0 {
        line_no += 1;

        if line_no % 1_000_000 == 0 {
            info!(
                "Processed {} lines",
                line_no.to_formatted_string(&Locale::en)
            );
        }

        if line.starts_with(b"#") {
            writer.write_all(&line)?;
            line.clear();
            continue;
        }

        let line_end = if line.ends_with(b"\n") {
            line.len() - 1
        } else {
            line.len()
        };
        let body_end = if line_end > 0 && line[line_end - 1] == b'\r' {
            line_end - 1
        } else {
            line_end
        };
        let body = &line[..body_end];

        let mut tab_count = 0usize;
        let mut attr_start = None;
        for (idx, &byte) in body.iter().enumerate() {
            if byte == b'\t' {
                tab_count += 1;
                if tab_count == 8 {
                    attr_start = Some(idx + 1);
                    break;
                }
            }
        }

        let Some(attr_start) = attr_start else {
            writer.write_all(&line)?;
            line.clear();
            continue;
        };

        // write column 1-8
        writer.write_all(&body[..attr_start])?;
        let mut wrote_attr = false;
        for raw_attr in body[attr_start..].split(|&b| b == b';') {
            let start = raw_attr
                .iter()
                .position(|b| !b.is_ascii_whitespace())
                .unwrap_or(raw_attr.len());
            let end = raw_attr
                .iter()
                .rposition(|b| !b.is_ascii_whitespace())
                .map(|idx| idx + 1)
                .unwrap_or(start);
            let attr = &raw_attr[start..end];
            if attr.is_empty() {
                continue;
            }

            let key_end = attr
                .iter()
                .position(|&b| b == b' ' || b == b'\t' || b == b'=')
                .unwrap_or(attr.len());
            let key = &attr[..key_end];
            let keep = keep_attrs.contains(&normalize_attr_key(key, args.keep_check_case));
            let is_isom_attr = if args.keep_check_case {
                key.starts_with(b"ISOM_")
            } else {
                key.len() >= 5 && key[..5].eq_ignore_ascii_case(b"ISOM_")
            };
            let remove = match args.chop_mode {
                ChopMode::All => !keep,
                ChopMode::Isomatch => is_isom_attr && !keep,
            };

            if remove {
                continue;
            }

            if wrote_attr {
                writer.write_all(b" ")?;
            }
            writer.write_all(attr)?;
            writer.write_all(b";")?;
            wrote_attr = true;
        }

        if !wrote_attr {
            writer.write_all(b".")?;
        }
        writer.write_all(&line[body_end..])?;
        line.clear();
    }

    writer.flush()?;

    info!("Chopped GTF saved to: {}", out_path.display());
    info!("Finished!");

    Ok(())
}

fn parse_keep_attrs(keep_attrs: Option<&str>, check_case: bool) -> HashSet<Vec<u8>> {
    let mut attrs: HashSet<Vec<u8>> = ["gene_id", "transcript_id"]
        .into_iter()
        .map(|attr| normalize_attr_key(attr.as_bytes(), check_case))
        .collect();

    attrs.extend(
        keep_attrs
            .into_iter()
            .flat_map(|attrs| attrs.split(','))
            .map(str::trim)
            .filter(|attr| !attr.is_empty())
            .map(|attr| normalize_attr_key(attr.as_bytes(), check_case)),
    );

    attrs
}

fn normalize_attr_key(key: &[u8], check_case: bool) -> Vec<u8> {
    if check_case {
        key.to_vec()
    } else {
        key.to_ascii_lowercase()
    }
}
