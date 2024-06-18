rm *.xml *.h5 *.out *.mp4
rm -rf __pycache__/
python *build.py
openmc .
python *plot.py
open *.mp4
