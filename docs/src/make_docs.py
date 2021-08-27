import os
import inspect

## modules to document
import firestudio.studios.gas_studio as gas_studio
import firestudio.studios.star_studio as star_studio

def write_doc(fname,module):

    fname = os.path.join(
        '../wiki',
        fname)

    with open(fname,'w') as handle:
        handle.write(module.__doc__)
    
def main():
    ## 
    write_doc('gas_studio.md',gas_studio)

    ##
    write_doc('star_studio.md',star_studio)
    

if __name__ == '__main__':
    main()
