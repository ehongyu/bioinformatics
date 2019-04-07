/*
JAVA APPLET TO GENERATE THE PHYLOGENETIC TREE OF A PROTEIN SEQUENCE FAMILY

Description:
Basically, it reads a multiple sequence alignment in CLUSTAL format,
and draws a phylogenetic tree graph in a web browser using UPGMA
algorithm.

Example:
To look at an example showing the working of this program, please go
to the following web link:
http://hongyu.org/software/upgma_example.html

Download:
Place to download the source code:
http://hongyu.org/software/upgma.tar.gz

Copyright:
Dr. Hongyu Zhang

Version:
1.1, Aug 17 2005

*/


import java.lang.*;
import java.applet.*;
import java.util.*;
import java.io.*;
import java.awt.*;
import java.net.URL;


public class Tree extends Applet {
    Node root;
    Node[] branches;
    Node[] finalBranches;
    float[][] distance; // distance matrix
    HashMap seqMap = new HashMap();
    int branchSize;    // initial branch size
    int newBranchSize;    // the present branch size
    String[] multipleSeqName;
    String[] multipleSeq;

    Image offscreenImg;
    Graphics offscreenG;

    public void init() {

        int i, seqNo;

        offscreenImg = createImage(this.getSize().width, this.getSize().height);
        offscreenG = offscreenImg.getGraphics();

        try{

            // read the alignment file in CLUSTALW format
            String alnFile = this.getParameter("clustalFile");
            BufferedReader in = new BufferedReader(new InputStreamReader(
                    new URL(alnFile).openStream()));

            //ignore the first three lines in the CLUSTALW output
            String s;
            for(i=0; i<3; i++) {
                s = in.readLine();
            }

            while((s=in.readLine()) != null) {
                if(s.length()<1 || s.charAt(0) == ' ') continue;

                // the seq line supposedly should have two columns,
                // one is the sequence name,
                // while the other is the amino acid or nucleotide
                String[] ss = s.split("\\s++");
                if(ss.length != 2) continue;

                if(seqMap.containsKey(ss[0])) {
                    seqMap.put(ss[0], (String)seqMap.get(ss[0]) + ss[1]);
                } else {
                    seqMap.put(ss[0], ss[1]);
                }
            }

            Object[] o = seqMap.keySet().toArray();
            Object[] o2 = seqMap.values().toArray();
            seqNo = o.length;

            multipleSeqName = new String[seqNo];
            multipleSeq = new String[seqNo];
            for ( i=0; i<seqNo; i++ ) {
                multipleSeqName[i] = (String)o[i];
                multipleSeq[i] = (String)o2[i];
            }

            // initialize branches
            branchSize = newBranchSize = seqNo;
            branches = new Node[branchSize];
            for(i=0; i<seqNo; i++) branches[i] = new Node(multipleSeqName[i], multipleSeq[i]);
            distance = new float[seqNo][seqNo];

            upgma();

        } catch (IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
        } catch (ArrayIndexOutOfBoundsException e) {
            System.err.println("Caught ArrayIndexOutOfBoundsException: " +
                    e.getMessage());
        }

    }


    //  generate tree

    public void upgma() {

        int i,j, minI, minJ;

        calculateDistance();

        while(newBranchSize > 1) { // repeat this cycle until newBranchSize=1
            float minDist = distance[0][1];
            minI = 0;
            minJ = 1;

            for(i=0; i<newBranchSize-1; i++) {
                for(j=i+1; j<newBranchSize; j++) {
                    if(distance[i][j] < minDist) {
                        minI = i;
                        minJ = j;
                        minDist = distance[i][j];
                    }
                }
            }

            refreshBranch(minI, minJ);
        }

        root = branches[0];

        // paint tree in the buffer

        int xUnit, yUnit, xShift, ticHeight;

        Font font = new Font("Helvetica", Font.BOLD, 10);
        offscreenG.setFont(font);

        setBackground(Color.white);
        setForeground(Color.black);

        // reorder the branches
        for(i=0; i<newBranchSize; i++) {

            Node tmp = branches[i];
            if(tmp.children.length != 0) {
                for(j=newBranchSize-1; j>i; j--) branches[j+1] = branches[j];
                branches[i] = tmp.children[0];
                branches[i+1] = tmp.children[1];
                newBranchSize++;
                i--;
            }
        }

      	for(i=0; i<branches.length; i++) multipleSeq[i] = branches[i].Sequence;

        calculateDistance();

        xUnit = 10;
        yUnit = 13;
        xShift = 5;
        ticHeight = 10;

        // draw background
        offscreenG.setColor(getBackground());
        offscreenG.fillRect(0, 0, this.getSize().width, this.getSize().height);
        offscreenG.setColor(getForeground());

        // set the position of leaves and draw the sequence name
        for(i=0; i<branchSize; i++) {
            offscreenG.drawString(branches[i].name, xUnit*50 + xShift*2, (int)(yUnit*(i+1.3)));
            branches[i].setPosition(xUnit*50 + xShift, yUnit*(i+1));
        }

        // draw axis and ticMarks
        offscreenG.drawLine(xShift, yUnit*(branchSize+1), xUnit*50 + xShift, yUnit*(branchSize+1));
        for(i=5; i<=10; i++) {
            offscreenG.drawLine(xUnit*10*(i-5) + xShift, yUnit*(branchSize+1),
                    xUnit*10*(i-5) + xShift, yUnit*(branchSize+1)+ticHeight);
            offscreenG.drawString("" + i*10, xUnit*10*(i-5),
                    yUnit*(branchSize+1) + 2*ticHeight);
        }

        // draw new nodes to branches. This is a similar cycle to the above
        // procedure.
        while(newBranchSize > 1) {

            float minDist = distance[0][1];
            minI = 0;
            minJ = 1;

            for(i=0; i<newBranchSize-1; i++) {
                for(j=i+1; j<newBranchSize; j++) {
                    if(distance[i][j] < minDist) {
                        minI = i;
                        minJ = j;
                        minDist = distance[i][j];
                    }
                }
            }

            refreshBranch(minI, minJ);

            float tmpx = xUnit * (50 - branches[minI].branchLength) + xShift;
            float tmpy = (branches[minI].children[0].y + branches[minI].children[1].y) / 2;
            branches[minI].setPosition( (int)tmpx, (int)tmpy );
            offscreenG.drawLine(branches[minI].children[0].x, branches[minI].children[0].y, branches[minI].x, branches[minI].children[0].y);
            offscreenG.drawLine(branches[minI].children[1].x, branches[minI].children[1].y, branches[minI].x, branches[minI].children[1].y);
            offscreenG.drawLine(branches[minI].x, branches[minI].children[0].y, branches[minI].x, branches[minI].children[1].y);

        }

        offscreenG.drawLine(branches[0].x, branches[0].y, xShift, branches[0].y);


    }

    public void paint(Graphics g) {
        System.out.println("painting ");
        g.drawImage(offscreenImg,0,0,this);
    }

    //    calculate distance matrix between branches[i]

    public void calculateDistance() {

        int seqNo = multipleSeq.length;
        int i, j, k, lengthAll;

        lengthAll = multipleSeq[0].length(); // the length of the alignment

        for(i=0; i<seqNo; i++) {
            for(j=0; j<seqNo; j++) {
                int length=0;
                int identity=0;
                for(k=0; k<lengthAll; k++) {
                    if(multipleSeq[i].charAt(k) == '-'
                            || multipleSeq[j].charAt(k) == '-') continue;
                    length++;
                    if(multipleSeq[i].charAt(k) == multipleSeq[j].charAt(k))
                        identity++;
                }

                distance[i][j] = 100 - (float)identity/length*100;
            }
        }

    }

    // generate new branch ii+jj, and remove the old branches ii & jj
    public void refreshBranch(int ii, int jj) {
        int i,j;

        // merge two children nodes to one parent node
        Node tmpNode = new Node(branches[ii], branches[jj]);
        float tmpf = distance[ii][jj] / 2;
        tmpNode.setBranchLength (tmpf);

        // replace branches[ii] with the new node
        branches[ii] = tmpNode;

        for(j=0; j<newBranchSize; j++) {
            distance[ii][j] = (distance[ii][j] + distance[jj][j]) / 2;
            distance[j][ii] = distance[ii][j];
        }

        // delete branches[jj]
        newBranchSize--;
        for(i=jj; i<newBranchSize; i++) {
            branches[i] = branches[i+1];
            for(j=0; j<newBranchSize+1; j++) distance[i][j] = distance[i+1][j];
        }
        for(j=jj; j<newBranchSize; j++) {
            for(i=0; i<newBranchSize; i++) distance[i][j] = distance[i][j+1];
        }

    }

    public void update() {
        System.out.println("updating..");
    }
}
