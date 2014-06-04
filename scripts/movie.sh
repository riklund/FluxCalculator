#!/bin/bash
if [ $# -ne 1 ]
then
	echo "Usage: $0 hist_dir"
fi

movieName="movie2.mkv"

ffmpeg -r 5 -i $1/%05d.png -c:v libx264 -preset veryslow -qp 0 -r 5 $movieName
