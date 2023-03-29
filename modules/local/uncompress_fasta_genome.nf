process UNCOMPRESS_FASTA_GENOME {
    tag "uncompress_fasta_genome"
    label 'process_high'

    input:
      path(f_ext)
    output:
      file "${base_name_file}.fasta"

    shell:
    //Tests to see if the input files are compressed. 
    //At this stage, only accepts fasta or .fa files, or .gz or .zip genome files.
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()
    if (ext_file == "fasta" | ext_file == "fa"){
    	'''
    	cp -n !{f_ext} !{base_name_file}.fasta
    	'''
    }else if(ext_file == "zip"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
  	  '''
  	  gunzip -f -S .zip !{f_ext}
  	  cp -n !{old_base_name_file} !{base_name_file}.fasta
  	  '''
    }else if(ext_file == "gz"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
	    gunzip -f !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else {
      '''
      echo "Your host genome files appear to have the wrong extension. \n Currently, the pipeline only supports .fasta or .fa, or compressed files with .zip or .gz extensions."
      '''
    }
}



