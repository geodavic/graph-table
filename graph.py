import networkx as nx
import numpy as np
from astropy.io import fits
import csv
import math

# TODO: arrow in middle of self-loops, alignment axis

class Polynomial(list): 
    # For printing polynomials -- credit to Eric on stack overflow
    def __repr__(self):
        # joiner[first, negative] = str
        joiner = {
            (True, True): '-',
            (True, False): '',
            (False, True): ' - ',
            (False, False): ' + '
        }

        result = []
        for power, coeff in reversed(list(enumerate(self))):
            j = joiner[not result, eval(coeff) < 0]

            coeff_str = str(abs(eval(coeff)))
            if "/" in coeff:
                [n,d] = coeff.split("/")
                n = int(n)
                d = int(d)
                coeff_str = f"\\frac{{ {abs(n)} }}{{ {abs(d)} }}"
            if eval(coeff) == 1 and power != 0:
                coeff_str = ''

            f = {0: '{}{}', 1: '{}{}x'}.get(power, '{}{}x^{}')

            result.append(f.format(j, coeff_str, power))

        return ''.join(result) or '0'


class TabularGraph:
    def __init__(self):
        self.matrices = self.get_matrices()
        self.subiso = self.get_lists("subiso.csv",cast_to_int=True)
        self.polys = self.get_lists("polys.csv")

    def get_matrices(self):
        hdul = fits.open("matrices.fits")
        mats = [np.array(hdul[i].data,np.int8) for i in range(len(hdul))]
        return mats

    def get_lists(self,file_name,cast_to_int=False):
        fp = open(file_name)
        reader = csv.reader(fp)
        lists = []
        for row in reader:
            if cast_to_int:
                lists.append([int(e) for e in row]) 
            else:
                lists.append(row)
        fp.close()
        return lists

    def _compute_loop_angle(self,px,py,loop_pos):
        angle_between = lambda v1,v2: np.arccos(np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))) 
        angle_from_e1 = lambda v1: angle_between(v1,np.array([1,0]))
        
        rel_pos = loop_pos - np.array([px,py])
        if rel_pos[1] < 0:
            angle = 2*math.pi - angle_from_e1(rel_pos)
        else:
            angle = angle_from_e1(rel_pos)

        angle = 360*(angle/(2*math.pi))
        return angle

    def _draw_line(self,px_out,py_out,px_in,py_in,bidir=False):
        middle_arrow_prop = 0.625
        line_str = f"\\draw[decoration={{markings, mark=at position {middle_arrow_prop} with {{\\arrow{{>}} }} }},"\
            +f"postaction={{decorate}}] ({px_out},{py_out}) -> ({px_in},{py_in});\n"
        if bidir:
            line_str += self._draw_line(px_in,py_in,px_out,py_out,bidir=False)
            line_str = line_str.replace("->","to[bend left]")
            
        return line_str

    def _draw_selfloop(self,px,py,loop_pos):
        distance = self.scale*9
        angle_width = 40
        opt_angle = self._compute_loop_angle(px,py,loop_pos)
        in_ = opt_angle + angle_width
        out_ = opt_angle - angle_width
        line_str = f"\\draw[->] ({px},{py}) edge[out={out_},in={in_},distance={distance}mm] ({px},{py});"
        return line_str 

    def _draw_subiso_poly(self,sub,poly,mid_y=0):
        if not sub:
            sub_str = "--"
        else:
            sub_str = str(sub)
        if not poly:
            poly_str = "--"
        else:
            poly_str = Polynomial(poly)
        return sub_str,poly_str

    def _pca_rotation(self,points):
        cov = np.cov(np.array(points).T)
        ev,eigv = np.linalg.eig(cov)

        max_idx = np.argmax(ev)
        [x,y] = eigv[:,max_idx]
        rot_mat = 0.7*self.scale*np.array([[x+y,y-x],[x-y,x+y]])
        return rot_mat

    def _get_layout(self,H,num_nodes,rounding=5):
        layout = nx.spring_layout(H,center=[0,0],seed=self.seed)
      
        # translate to centroid of non-self-loop nodes
        centroid = sum([v for k,v in layout.items() if k<num_nodes])/num_nodes
        layout = {k:v-centroid for k,v in layout.items()}

        # scale so it fits in unit circle
        sc = max([np.linalg.norm(v) for k,v in layout.items() if k<num_nodes])
        layout = {k:v/sc for k,v in layout.items()}

        # rotate to principal axes (this rotates to 45deg)
        rot = self._pca_rotation(list(layout.values()))

        # set alignment (if something other than 45 wanted)
        cos = math.cos(self.align_angle)
        sin = math.sin(self.align_angle)
        rot_align = np.array([[cos,-sin],[sin,cos]])

        trans = np.matmul(rot_align,rot)
        layout = {k:trans@v for k,v in layout.items()}

        layout = {k:np.round(v,decimals=rounding) for k,v in layout.items()}
        return layout

    def create_block(self,mat,sub,poly):
        hshift = 0
        vshift = 0
        H = nx.from_numpy_matrix(mat)
        G = nx.from_numpy_matrix(mat)

        # hanlde self-loops
        for idx in np.nonzero(np.diag(mat))[0]:
            H.remove_edge(idx,idx) 
            H.add_edge(idx,len(mat)+idx)

        layout = self._get_layout(H,len(mat))

        block_str = ""
        x_vals = []
        y_vals = []

        # draw nodes of graph
        for node in range(len(mat)):
            pos = layout[node]
            px = pos[0] + hshift
            py = pos[1] + vshift
            x_vals.append(px)
            y_vals.append(py)
            block_str += f"\\node[shape=circle,draw=black,fill=black,scale=0.3] at ({px},{py}) {{}};\n"

        # draw ghost nodes for bounding figure
        for node_x in [-self.scale,self.scale]:
            for node_y in [-self.scale,self.scale]:
                block_str += f"\\node at ({node_x},{node_y}) {{}};\n";

        asym_idx = np.where((mat - mat.T)>0) # unidirectional edges
        sym_idx = np.nonzero(np.triu(mat*mat.T)) # bidirectional edges (or self-loops)

        # draw unidirectional edges
        for edge_out,edge_in in zip(*asym_idx):
            px_out = layout[edge_out][0] + hshift
            py_out = layout[edge_out][1] + vshift
            px_in = layout[edge_in][0] + hshift
            py_in = layout[edge_in][1] + vshift
            block_str += self._draw_line(px_out,py_out,px_in,py_in,bidir=False)

        # draw bidirectional edges
        for edge_out,edge_in in zip(*sym_idx):
            px_out = layout[edge_out][0] + hshift
            py_out = layout[edge_out][1] + vshift
            px_in = layout[edge_in][0] + hshift
            py_in = layout[edge_in][1] + vshift
            if edge_out != edge_in:
                block_str += self._draw_line(px_out,py_out,px_in,py_in,bidir=True)
            else:
                # draw self-loops
                loop_pos = layout[len(mat)+edge_out]
                block_str += self._draw_selfloop(px_in,py_in,loop_pos)

        block_str = "\\multicolumn{1}{m{2.5cm}}{\\begin{tikzpicture}\n" + block_str + "\\end{tikzpicture}}\n"
        return block_str

    def make_tabular(self,tikz_str):
        titles = " \small{Graph} & \small{Subisomorphisms} & \small{Polynomial} "
        top = "\\hline \n"+titles+"&"+titles+"\\\\ \n \\hline \n"
        rval = "\\begin{tabular}{ccc||ccc}\n"+ top + tikz_str + "\\end{tabular}\n"
        return rval

    def render(self,num=None,seed=1,scale=1,angle=45,out_name="input.tex"):
        self.seed = seed
        self.scale = scale
        self.align_angle = 2*math.pi*((angle-45)/360)

        tikz_str = ""
        counter = True
        if not num:
            num = len(self.matrices)
        for mat,sub,poly in list(zip(self.matrices,self.subiso,self.polys))[:num]:
            block_str = self.create_block(mat,sub,poly)
            sub_str,poly_str = self._draw_subiso_poly(sub,poly)
            tikz_str += block_str
            if counter:
                tikz_str += f"& ${sub_str}$ & ${poly_str}$ &\n"
            else:
                tikz_str += f"& ${sub_str}$ & ${poly_str}$ \\\\ \n"
            counter ^= 1
       
        final = self.make_tabular(tikz_str)

        with open(out_name,"w") as fp:
            fp.write(final)
        print(f"Output written to {out_name}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser() 
    parser.add_argument("--scale",type=float,default=1,help="Scale factor for graphs (bigger value -> bigger graphs)")
    parser.add_argument("--seed",type=int,default=1,help="Seed for graph drawing function.")
    parser.add_argument("--limit",type=int,default=None,help="Max number of graphs to render.")
    parser.add_argument("--angle",type=float,default=45,help="Alignment angle of graphs (in degrees).")

    args = parser.parse_args()
    T = TabularGraph()
    T.render(scale=args.scale,seed=args.seed,num=args.limit,angle=args.angle)