# ffmpeg -framerate 30 -pattern_type glob -i 'fig/*.png' \
#   -c:v libx264 -pix_fmt yuv420p 'fig/out.mp4'

magick -delay 20 -loop 0 'fig/*.png' 'fig/out.gif'