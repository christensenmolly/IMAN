import PIL
from PIL import Image

import sys



def main(input_image, output_image, dpi=100):
    basewidth = dpi
    img = Image.open(input_image)

    wpercent = (basewidth / float(img.size[0]))
    hsize = int((float(img.size[1]) * float(wpercent)))
    img = img.resize((basewidth, hsize), PIL.Image.ANTIALIAS)
    img.save(output_image)



if __name__ == "__main__":
     input_image = sys.argv[1]
     output_image = sys.argv[2]
     dpi = sys.argv[3]
     main(input_image, output_image, dpi=int(dpi))