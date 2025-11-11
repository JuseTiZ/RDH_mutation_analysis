mkdir -p utils figures outputs
for fig in {1..5}; do
    mkdir -p figures/Fig$fig
done

for fig in Fig*/figure*; do

    newname=$(echo $fig | sed 's/figure/fig/')
    mv $fig $newname

done

git init
git add *
git commit -m "init commit"
git branch -M main
git remote add origin https://github.com/JuseTiZ/RDHGs_mutation_analysis.git
git push -u origin main
