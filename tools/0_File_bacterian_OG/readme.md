
This step should only be performed by the programmer each time OrthoDB is updated. It allows only bacterial ortholog groups to be extracted from the OrthoDB database.
```bash=
python formatting_bacterial_orthologue_file_final.py -o $PATH_ORTHODB/Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
```
