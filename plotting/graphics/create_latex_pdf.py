import matplotlib
import sys
from pylatex import Document, Package, Section, Figure, NoEscape, SubFigure
from pylatex.utils import italic, bold


matplotlib.use('Agg')  # Not to use X server. For TravisCI.
import matplotlib.pyplot as plt  # noqa


def main(name, pictures):
    doc = Document()
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))

    doc.append(bold(name))
    




    '''
    with doc.create(Figure(position='h!')) as pic:
            pic.add_image(pictures[0], width=NoEscape(r'0.45\linewidth'))
            pic.add_caption('Picture of %s' % (name))

    with doc.create(Figure(position='h!')) as kittens:
            with doc.create(SubFigure(
                    position='c',
                    width=NoEscape(r'0.45\linewidth'))) as left_kitten:

                left_kitten.add_image(pictures[1],
                                      width=NoEscape(r'\linewidth'))
                left_kitten.add_caption('IRAF results')

    with doc.create(Figure(position='h!')) as kittens:
            with doc.create(SubFigure(
                    position='c',
                    width=NoEscape(r'0.45\linewidth'))) as right_kitten:

                right_kitten.add_image(pictures[2],
                                       width=NoEscape(r'\linewidth'))
                right_kitten.add_caption('Horizontal profile')
            with doc.create(SubFigure(
                    position='c',
                    width=NoEscape(r'0.45\linewidth'))) as right_kitten:

                right_kitten.add_image(pictures[3],
                                       width=NoEscape(r'\linewidth'))
                right_kitten.add_caption('Vertical profile')
            kittens.add_caption("Profiles")

    with doc.create(Figure(position='h!')) as kittens:
            with doc.create(SubFigure(
                    position='c',
                    width=NoEscape(r'0.45\linewidth'))) as right_kitten:

                right_kitten.add_image(pictures[4],
                                       width=NoEscape(r'\linewidth'))
                right_kitten.add_caption('Model')


    '''

    doc.generate_pdf('test')
    

if __name__ == '__main__':
    name = str(sys.argv[1]) 
    pictures = str(sys.argv[2])
    
    try:
        pictures = pictures.split(',')
    except:
        pictures = [pictures]
    main(name, pictures)