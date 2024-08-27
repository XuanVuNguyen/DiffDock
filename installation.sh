conda install \
    pytorch==1.13.0 \
    pytorch-cuda=11.7 \
    -c pytorch -c nvidia
pip install \
    torch-scatter==2.0.9 \
    torch-sparse==0.6.15 \
    torch-cluster==1.6.0 \
    torch-spline-conv==1.2.1 \
    torch-geometric==2.0.4 \
    rdkit==2022.03.3 \
    e3nn==0.5.0 \
    fair-esm==2.0.0 \
    biopython==1.84 \
    -f https://data.pyg.org/whl/torch-1.13.0+cu117.html