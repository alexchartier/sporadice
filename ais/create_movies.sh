#!/bin/zsh 
# from directory inside plot directory

a=1
mkdir tmp
outdir=~/data/sporadice/plots/ais_era5_movie.mp4
for i in  ~/data/sporadice/plots/ais_era5_maps/*png; do 
    new=$(printf "tmp/%04d.png" "$a")
    cp "$i" "$new"
    let a=a+1  
done

ffmpeg  -f image2 -r 6 -pattern_type glob -i 'tmp/*.png' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"  -pix_fmt yuv420p $outdir
rm -rf tmp


