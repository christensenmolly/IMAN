class Settings():
    def __init__(self):

        self.known_catalogs = ['SDSS', 'PS1', 'PS', 'NOMAD', 'NOMAD1']



    def _info(self):
        s_dict = self.__dict__
        for name in s_dict:
            print('%s: %s' %(name, s_dict[name]))

if __name__ == '__main__':
    s = Settings()
    s._info()
        
