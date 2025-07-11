# ffmpeg -framerate 30 -pattern_type glob -i 'scripts/fig/*.png' \
#   -c:v libx264 -pix_fmt yuv420p scripts/fig/out.mp4

magick -delay 20 -loop 0 'scripts/fig/*.png' 'scripts/fig/out.gif'