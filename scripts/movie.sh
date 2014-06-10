#!/bin/bash
if [ $# -ne 2 ]
then
	echo "Usage: $0 hist_dir movie_name.mkv"
	exit 1;
fi

movieName=$2

ffmpeg -r 5 -i $1/%05d.png -c:v libx264 -preset veryslow -qp 0 -r 5 $movieName
