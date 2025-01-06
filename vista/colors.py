import matplotlib.pyplot as plt

class colors:

    '''
    Colors class:reset all colors with colors.reset; two
    sub classes fg for foreground
    and bg for background; use as colors.subclass.colorname.
    i.e. colors.fg.red or colors.bg.greenalso, the generic bold, disable,
    underline, reverse, strike through,
    and invisible work with the main class i.e. colors.bold
    '''
    reset         = '\033[0m'
    bold          = '\033[01m'
    disable       = '\033[02m'
    underline     = '\033[04m'
    reverse       = '\033[07m'
    strikethrough = '\033[09m'
    invisible     = '\033[08m'

    class fg:
        black      = '\033[30m'
        red        = '\033[31m'
        green      = '\033[32m'
        orange     = '\033[33m'
        blue       = '\033[34m'
        purple     = '\033[35m'
        cyan       = '\033[36m'
        lightgrey  = '\033[37m'
        darkgrey   = '\033[90m'
        lightred   = '\033[91m'
        lightgreen = '\033[92m'
        yellow     = '\033[93m'
        lightblue  = '\033[94m'
        pink       = '\033[95m'
        lightcyan  = '\033[96m'

    class bg:
        black      = '\033[40m'
        red        = '\033[41m'
        green      = '\033[42m'
        orange     = '\033[43m'
        blue       = '\033[44m'
        purple     = '\033[45m'
        cyan       = '\033[46m'
        lightgrey  = '\033[47m'
    

# ----------------------------------------------------------------------
# >>> Palette                                    
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/01/06  - created
#
# Desc:
#     customized colors for scientific plotting
# ----------------------------------------------------------------------

class Palette:
    
    """
    A class to store and manage custom colors for easy access.
    """
    
    def __init__(self):
        
        # Define custom colors as a dictionary
        
        self.colors = {
            "blue"        : "#1f77b4",
            "orange"      : "#ff7f0e",
            "green"       : "#2ca02c",
            "red"         : "#d62728",
            "purple"      : "#9467bd",
            "cyan"        : "#17becf",
            "yellow"      : "#bcbd22",
            "dark_gray"   : "#7f7f7f",
            "steelblue"   : "#2878b5",
            "skyblue"     : "#9ac9db",
            "lightsalmon" : "#f8ac8c",
            "brickred"    : "#c82423",
            "lightcoral"  : "#ff8884",
        }

    def add_color(self, name, hex_code):
        """Add a custom color to the dictionary."""
        if name in self.colors:
            raise ValueError(f"Color '{name}' already exists.")
        self.colors[name] = hex_code

    def get(self, name):
        """Retrieve a color by name."""
        if name not in self.colors:
            raise KeyError(f"Color '{name}' not found.")
        return self.colors[name]

    def show_colors(self):
        """show all available colors."""

        n_colors = len( self.colors )
        
        fig, ax = plt.subplots( figsize=(8, n_colors*0.5) )
        
        for i, (name, hex_code) in enumerate(reversed(self.colors.items())):
            ax.add_patch(plt.Rectangle((0,i),1,1,color=hex_code))
            ax.text(1.1,i+0.5,f"{name} ({hex_code})",va='center',fontsize=10)
        
        plt.bar([5.5,55,550],[1.0,2.0,3.0],color=['r','b','g'],alpha=0.5)
        ax.set_xlim(0,2)
        ax.set_ylim(0,n_colors)
        ax.axis("off")
        plt.title("Palette Choices", fontsize=14, weight='bold')
        plt.tight_layout()
        plt.show()


# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2025/01/06  - created
#
# Desc
#
# ----------------------------------------------------------------------


if __name__ == "__main__":

    palette = Palette()
    palette.show_colors()

        