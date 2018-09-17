#!/bin/bash
CUR_DIR="./scrap_results/sou/"
for x in ./presidential_speech/*; do
	if ! [ -d "$x" ]; then
		if [ "$x" == "./presidential_speech/scrapy.cfg" ]; then
			echo "is scrap"
		else
			mv -- "$x" "$CUR_DIR"
		fi
	fi
done
