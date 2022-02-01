#!/bin/bash

xclip -selection clipboard -target image/png -out > Screenshot_"$(date "+%Y-%m-%d_%H-%M-%S")".png
