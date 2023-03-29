process UNCOMPRESS_GFF {
    tag "uncompress_GFF"
    label 'process_high'

    input:
        path(f_ext) 
    output:
        file "${base_name_file}.gff3"

    shell:
    //Tests to see if the input files are compressed. 
    //At this stage, only accepts .gff or gff3, or .gz or .zip annotation files.
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()

    if (ext_file == "gff" | ext_file == "gff3"){
      '''
      cp -n !{f_ext} !{base_name_file}.gff3
      '''
    }else if(ext_file == "zip"){
      '''
      gunzip -f -S .zip !{f_ext}
      cp -n !{base_name_file} !{base_name_file}.gff3
      '''
    }else if(ext_file == "gz"){
      //if gff or gff3, need to save as .gff3
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.gff|.gff3/,"")
      '''
      gunzip -f !{f_ext}
      cp -n !{old_base_name_file} !{base_name_file}.gff3
      '''
    }else {
      '''
      echo "Your GFF file appears to be in the wrong format or has the wrong extension. \n Currently, the pipeline only supports .gff or .gff3, or compressed files with .zip or .gz extensions."
      '''
    }
}
