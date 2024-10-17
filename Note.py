
# variable list for INCA's output with second order correlation.

'''
VarList  = ['x',                'y',                    '<u>',                  \
            '<v>',              '<w>',                  '<rho>',                \
            '<rho*E>',          '<p>',                  '<T>',                  \
            '<mu>',             '||grad(<velocity>)||', '||curl(<velocity>)||', \
            'div(<velocity>)',  '||shear(<velocity>)||','||grad(<T>)||',        \
            '<u`u`>',           '<u`v`>',               '<u`w`>',               \
            '<u`rho`>',         '<u`rho*E`>',           '<u`p`>',               \
            '<u`T`>',           '<u`mu`>',              '<v`v`>',               \
            '<v`w`>',           '<v`rho`>',             '<v`rho*E`>',           \
            '<v`p`>',           '<v`T`>',               '<v`mu`>',              \
            '<w`w`>',           '<w`rho`>',             '<w`rho*E`>',           \
            '<w`p`>',           '<w`T`>',               '<w`mu`>',              \
            '<rho`rho`>',       '<rho`rho*E`>',         '<rho`p`>',             \
            '<rho`T`>',         '<rho`mu`>',            '<rho*E`rho*E`>',       \
            '<rho*E`p`>',       '<rho*E`T`>',           '<rho*E`mu`>',          \
            '<p`p`>',           '<p`T`>',               '<p`mu`>',              \
            '<T`T`>',           '<T`mu`>',              '<mu`mu`>',             \
            'walldist',         'shocksens',            'entropy']
'''

# compile the fortran library in vista.lib

# >>> cd vista/lib/
# >>> f2py -m form -c *.f90