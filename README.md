# ROTMOL


PROGRAMA PARA ROTAÇÃO DE MOLÉCULAS EM TORNO DA ORIGEM

COMPILAÇÃO: gfortran -O2 -o rotmol rotmol.f90

INPUT:
3 3              => numero de átomos | numero de operações de rotação
10 5 7 3         => numero de repetições da 1a op | phi | theta | gamma
30 1 0 2         => numero de repetições da 2a op | phi | theta | gamma
90 2 2 2         => numero de repetições da 3a op | phi | theta | gamma
                 => LINHA EM BRANCO
H 0.0 0.0 1.0    => COORDENADAS DOS N ÁTOMOS
O 0.0 0.0 0.0
H 0.0 1.0 0.0


A cada REPETIÇÃO será salvo a estrutura rotacioanda... no exemplo acima serão gerados 130 imagens + 1 (inicial) = 131 imagens.

INFO:  DEBUG => colocar o número de átomos negativo (-3, por ex.)


//TODO: Melhorar a escrita do DEBUG
