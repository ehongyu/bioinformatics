/*
  JAVA APPLET TO GENERATE THE PHYLOGENETIC TREE OF A PROTEIN SEQUENCE FAMILY

*/


import java.lang.*;

class Node {
    String name;
    String Sequence;
    Node parent;
    Node[] children;
    float branchLength;
    int x, y;


    public Node() {
        name = "Noname";
    };

    public Node(String str1, String str2) {     //leaf
        name = str1;
        Sequence = str2;
        children = new Node[0];
        branchLength = 0;
    }

    public Node(Node n0, Node n1) {
        this.children = new Node[2];
        children[0] = n0;
        children[1] = n1;
        n0.setParent(this);
        n1.setParent(this);
    }

    public void setChildren(Node n0, Node n1) {
        children = new Node[2];
        children[0] = n0;
        children[1] = n1;
        n0.setParent(this);
        n1.setParent(this);
    }

    public void setParent(Node np) {
        parent = np;
    }

    public void setBranchLength(float tmp) {
        branchLength = tmp;
    }

    public void setPosition(int i, int j) {
        x = i;
        y = j;
    }
}
