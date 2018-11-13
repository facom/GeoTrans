#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#IPYTHON SPECIFICS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#IPython
from IPython.display import HTML, Image
import IPython.core.autocall as autocall
class ExitClass(autocall.IPyAutocall):
    '''
         Simulation autcall clas
    '''
    def __call__(self):
        exit()

class dictObj(object):
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        self.__dict__.update(other.__dict__)
        return self
        
def in_ipynb():
    try:
        cfg = get_ipython().config 
        return True
    except NameError:
        return False

QIPY=False
if in_ipynb():
    QIPY=True
if not QIPY:
    def Image(url="",filename="",f=""):pass
    def get_ipython():
        foo=dictObj(dict())
        foo.run_line_magic=lambda x,y:x
        foo.magic=lambda x:x
        return foo

#GRAPHICAL
if not QIPY:
    from matplotlib import use
    use('Agg')
import matplotlib.pyplot as plt
