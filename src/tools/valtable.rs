use log::info;

use crate::{ValTableArgs, tools::tools_error::ToolError};




pub fn run_valtable(args:&ValTableArgs) -> Result<(), ToolError> {


    // load merged gtf check if it has ISOM headers, parse the file name and the sample id string(S1 etc) 

    // and build the projection of (src_tx_id,file_id) to (isomatch_tx_id, isomatch_tx index in the file) 

    // read the first file path, remove path, only get the file name get the file idx, if failed, then the 
    // skip this file and warning

    // read the attribute from transcript line and take the target attribute (args.attr_val) and transcript_id
    // make a Vec with default value (args.default_value), length is the number of isomatch_tx in the merged gtf file. 

    // look up the (isomatch_tx_id, isomatch_tx index in the file)  based on  (src_tx_id,file_id) 
    // fill up the Vec with attribute value

    // write it temporarly in to the disk 

    // after all files processed 

    // read all the temporary files and generate a matrix that have the first columns for isomatch transcirpt id 
    // and each file for one rest of the columns. 

    // you need make a stats struct to log the status of each file, how many transcirpts been successfuly exttacted,
    // how many failed. 

    info!("Finished!");
    Ok(())
}