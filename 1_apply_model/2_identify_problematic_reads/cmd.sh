# Get readnames with repeat counts greater than or equal to 10
awk -F "\t" '{if($2>=10) print}' chm13.nanopore.rel7.map_chm13.sorted.allchromosomes.repeatsTotal > chm13.nanopore.rel7.map_chm13.sorted.allchromosomes.repeatsTotal.gt10
cut -f1 chm13.nanopore.rel7.map_chm13.sorted.allchromosomes.repeatsTotal.gt10 > chm13.nanopore.rel7.map_chm13.sorted.allchromosomes.repeatsTotal.gt10.readnames
