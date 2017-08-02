README.md of /blobology/dev

––––––––––––––––––––––––––––––

Input formats

legacy format:

| columns |    ID   |  len |   gc  |  cov1 |  cov2 |  cov3 |     taxlevel1    |    taxlevel2   |    taxlevel3    |       taxlevel4       |
|:-------:|:-------:|:----:|:-----:|:-----:|:-----:|:-----:|:----------------:|:--------------:|:---------------:|:---------------------:|
|   type  |   str   |  int | float | float | float | float |        str       |       str      |       str       |          str          |
| example |    ID   |  len |   gc  |  cov1 |  cov2 |  cov3 | taxlevel_species | taxlevel_order | taxlevel_phylum | taxlevel_superkingdom |
|         | contig1 | 1000 |  0.56 |  100  | 150.2 |  80.3 |    Danio rerio   |  Cypriniformes |     Chordata    |       Eukaryota       |
|         | contig2 |  500 |  0.54 |  12.9 |  11.1 |  15.7 |   Not Annotated  |  Not Annotated |  Not Annotated  |     Not Annotated     |


- eval format as in gc_cov_annotate_v1.pl --evalue (requires BLAST -outfmt '6 qseqid taxid std')

| columns | ID      | len  | gc    | cov1  | cov2  | cov3  | taxlevel1        | taxlevel2      | taxlevel3       | taxlevel4             | eval       |
|---------|---------|------|-------|-------|-------|-------|------------------|----------------|-----------------|-----------------------|------------|
| type    | str     | int  | float | float | float | float | str              | str            | str             | str                   | scientific |
| example | ID      | len  | gc    | cov1  | cov2  | cov3  | taxlevel_species | taxlevel_order | taxlevel_phylum | taxlevel_superkingdom | eval       |
|         | contig1 | 1000 | 0.56  | 100   | 150.2 | 80.3  | Danio rerio      | Cypriniformes  | Chordata        | Eukaryota             | 1e-25      |
|         | contig2 | 500  | 0.54  | 12.9  | 11.1  | 15.7  | Not Annotated    | Not Annotated  | Not Annotated   | Not Annotated         | N/A        |

- novel format (to be inplemented)

| columns |    ID   |  len |   gc  |              cov              |                                            tax                                            | eval       |
|:-------:|:-------:|:----:|:-----:|:-----------------------------:|:-----------------------------------------------------------------------------------------:|------------|
|   type  |   str   |  int | float |             float             |                                            str                                            | scientific |
| example |    ID   |  len |   gc  |              cov              |                                            tax                                            | eval       |
|         | contig1 | 1000 |  0.56 | cov1=100;cov2=150.2;cov3=80.3 |       species=Danio rerio;order=Cypriniformes;phylum=Chordata;superkingdom=Eukaryota      | 1e-25      |
|         | contig2 |  500 |  0.54 | cov1=12.9;cov2=11.1;cov3=15.7 | species=Not Annotated;order=Not Annotated;phylum=Not Annotated;superkingdom=Not Annotated | N/A        |
