protein_ligand_csv=mchr1/diffdock_input.csv
out_dir=results_October/rosettafold_active/

cp -r $protein_ligand_csv $out_dir

python -m inference \
    --config default_inference_args.yaml  \
    --protein_ligand_csv $protein_ligand_csv \
    --out_dir $out_dir