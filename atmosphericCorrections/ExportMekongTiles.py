import os

pyPath = '/usr/bin/python'

cmd = pyPath +' ExportLandsatSRComposite.py \
        -y 2002 \
        -p ./eeAtsCorLut/ \
        -s drycool \
        -b 96 21 98 23'

os.system(cmd)
