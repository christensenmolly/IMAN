from PIL import Image
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

def main(input_png_1, input_png_2, output_image='merged.png', verbosity=True):
    if verbosity: print('Merging two images...')

    # Take two images for blending them together   
    image1 = Image.open(input_png_1)
    image2 = Image.open(input_png_2)


    ny = np.shape(image1)[0]
    nx = np.shape(image1)[1]

    # Make sure images got an alpha channel
    image5 = image1.convert("RGBA")
    image6 = image2.convert("RGBA")



    # alpha-blend the images with varying values of alpha
    #alphaBlended1 = Image.blend(image5, image6, alpha=.2)
    alphaBlended2 = Image.blend(image5, image6, alpha=.4)

    # Display the alpha-blended images
    #alphaBlended1.show()
    #alphaBlended2.show()
    if len(output_image.split('.png'))<2:
        alphaBlended2.save('tmp.png')
    else:
        alphaBlended2.save(output_image)
        return 0
    fig =  plt.figure(0)
    ax = plt.subplot()
    
    img=mpimg.imread('tmp.png')
    imgplot = plt.imshow(img)

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    

    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()
    os.remove('tmp.png')
    
    if verbosity: print('Done!')
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merging two rgb images")
    parser.add_argument("input_png1", help="First input png image", type=str)
    parser.add_argument("input_png2", help="Second input png image", type=str)
    parser.add_argument("--o", help="Output image", type=str, default='merged.png')    
    args = parser.parse_args()
    
    input_png1 = args.input_png1
    input_png2 = args.input_png2
    output_image = args.o
    main(input_png1, input_png2, output_image) 