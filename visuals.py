from tkinter import *


def vis_brak(tup_list, sequence):
    rtn = []
    for i in sequence:
        rtn.append(".")
    for b1, b2 in tup_list:
        rtn[b1] = "("
        rtn[b2] = ")"
    return "".join(rtn)


def vis_arc(tup_list, sequence):
    root = Tk()
    root.resizable(True, True)
    c_h = 960
    c_w = 1280
    c = Canvas(root, height=c_h, width=c_w, bg="white")
    c.pack()
    # create_line(x1,y1,x2,y2)
    c.create_line(100, c_h-100, c_w-100, c_h-100, width=2)
    n = len(sequence)
    space = (c_w-200)/(n-1)
    loc = {}
    seq_len = len(sequence)
    for i in range(seq_len):
        x = (space*i)+100
        y = c_h-100
        if seq_len < 200:
            c.create_text(x, y+15, text=sequence[i])
            c.create_line(x, y+5, x, y-5)
        loc[i] = (x, y)
    for b1, b2 in tup_list:
        # create_arc(coord bounding box(x1,y1 to x2,y2),start ang, extend ang)
        c.create_arc(loc[b1][0], loc[b1][1]-(b2-b1)*500*(1/len(sequence)),
                     loc[b2][0], loc[b2][1]+(b2-b1)*500*(1/len(sequence)),
                     start=0, extent=180, style="arc", width=2)

    root.mainloop()
