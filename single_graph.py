from graph import TabularGraph

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser() 
    parser.add_argument("--scale",type=float,default=1,help="Scale factor for graphs (bigger value -> bigger graphs)")
    parser.add_argument("--seed",type=int,default=1,help="Seed for graph drawing function.")
    parser.add_argument("--angle",type=float,default=45,help="Alignment angle of graphs (in degrees).")
    parser.add_argument("--self_loop_nodes",action="store_true",help="Treat self-loops as nodes when making each graph layout.")
    parser.add_argument("--file",type=str,help="The adjacency matrix file. Can be either `.fits` or `.csv`")

    args = parser.parse_args()
    T = TabularGraph()

    T.single_graph_render(args.file,scale=args.scale,seed=args.seed,angle=args.angle,loops_are_nodes=args.self_loop_nodes)
