mkdir -p utils figures outputs
for fig in {1..5}; do
    mkdir -p figures/Fig$fig
done

for fig in Fig*/figure*; do

    newname=$(echo $fig | sed 's/figure/fig/')
    mv $fig $newname

done

git init
git add -A
git commit -m "update"
git push