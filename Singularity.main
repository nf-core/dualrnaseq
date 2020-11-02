Bootstrap: docker
From: nfcore/base


%environment

    PATH=/opt/conda/envs/nf-core-dualrnaseq-1.0dev/bin:$PATH
    export PATH

    PYTHONPATH=/opt/conda/envs/nf-core-dualrnaseq-1.0dev/lib/python3.7/site-packages/:$PYTHONPATH
    export PYTHONPATH


%files
    environment.yml 

%post

    alias conda="/opt/conda/bin/conda"    

    /opt/conda/bin/conda env create -f /environment.yml
    conda clean --tarballs --index-cache --source-cache 
