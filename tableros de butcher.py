from pylab import *

gauss2 = array(['Gauss 2 etapas', array([[1/4, 1/4 - sqrt(3)/6], [1/4 + sqrt(3)/6, 1/4]]), array([1/2, 1/2])], dtype = 'O')
gauss3 = array(['Gauss 3 etapas', array([[5/36, 2/9 - sqrt(15)/15, 5/36 - sqrt(15)/30],[5/36 + sqrt(15)/24, 2/9, 5/36 - sqrt(15)/24],[5/36+sqrt(15)/30, 2/9 + sqrt(15)/15, 5/36]]), array([5/18, 4/9, 5/18])], dtype='O')

radau1a2 = array(['Radau IA 2 etapas', array([[1/4, -1/4], [1/4, 5/12]]), array([1/4, 3/4])], dtype = 'O')
radau2a2 = array(['Radau IIA 2 etapas', array([[5/12, -1/12], [3/4, 1/4]]), array([3/4, 1/4])], dtype = 'O')
radau2a3 = array(['Radau IIA 3 etapas', array([[(88 - 7*sqrt(6))/360, (296-169*sqrt(6))/1800, (-2+3*sqrt(6))/225], [(296+169*sqrt(6))/1800, (88+7*sqrt(6))/360, (-2-3*sqrt(6))/225],[(16-sqrt(6))/36,(16+sqrt(6))/36, 1/9]]), array([(16-sqrt(6))/36, (16+sqrt(6))/36,1/9])], dtype='O')

lobatto3c2 = array(['Lobatto IIIC 2 etapas', array([[1/2, -1/2], [1/2, 1/2]]), array([1/2, 1/2])], dtype = 'O')
lobatto3c3 = array(['Lobatto IIIC 3 etapas', array([[1/6,-1/3,1/6],[1/6,5/12,-1/12],[1/6,2/3,1/6]]), array([1/6,2/3,1/6])], dtype = 'O')


