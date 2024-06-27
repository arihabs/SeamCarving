import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;
import edu.princeton.cs.algs4.StdOut;
//import java.util.TreeSet;
import java.util.Arrays;
import java.util.Stack;
//TODO:
// Need to keep track if this.H and this.W change after transposing and look at methods that rely on these variables.
// Need to update energies of pixels that bordered pixels that were removed.
public class SeamCarver {
    private Picture pic;
    private double[][] pixEnergyArr;

    private double[][] pixEnergyArr2;

    private int[][] pixels;

    private boolean isTransposed = false;
    private boolean seamRemoved = false;
    private final double BORDER_ENERGY = 1000;
    private int W;
    private int H;

    private int nPixels(){
        return W*H;
    }

    private int nVertices(){
        return nPixels()+2;
    }

    private int srcVertex = 0;
    private int sinkVertex = nVertices()-1;

    private int[] edgeTo;
    private double[] distTo;

    // Keep set of existing pixels/vertices (ones that haven't been deleted).
//    private TreeSet<Integer> existingPixelSet;


    private int[] rgbIntToArray(int rgb){
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >>  8) & 0xFF;
        int b = (rgb >>  0) & 0xFF;
        int[] arr = {r,g,b};
        return arr;
    }

    private int idxOffset = 1;
    private int[] idxToColRow(int i){
        int col = i / W;
        int row = i % W;
//        int col = (i-idxOffset) / W;
//        int row = (i-idxOffset) % W;
        int[] colRow = {col,row};
        return colRow;
    }
    private int[] vertexToColRow(int i){
        return idxToColRow(i-idxOffset);
//        int col = (i-idxOffset) / W;
//        int row = (i-idxOffset) % W;
//        int[] colRow = {col,row};
//        return colRow;
    }
    private int colrowToIndex(int col, int row){
//        return row*W + col + idxOffset;
        return row*W + col;
    }

    private int colrowToVertex(int col, int row){
        return colrowToIndex(col, row) + idxOffset;
    }

    private double computePixEnergy(int col, int row){
        if(col == 0 || row == 0 || col == this.W-1 || row == this.H-1)
            return BORDER_ENERGY;

        double ex = energyCol(col, row,false);
        double ey = energyCol(col, row,true);
        return Math.sqrt(ex+ey);
    }
    private double energyCol(int col, int row, boolean transpose){
//        assert(col > 0 && col < width()-1);
//        assert(row > 0 && col < height()-1);
        double ex = 0.0;
//        int[] rgbArrCurr = rgbIntToArray(pic.getRGB(col,row));
//        int[] rgbArrCurr = rgbIntToArray(pixels[col][row]);
        int[] rgbArrCurr = rgbIntToArray(pixels[row][col]);
        int[] rgbArrPrev;
        int[] rgbArrNext;
        if(!transpose) {
//            rgbArrPrev = rgbIntToArray(pic.getRGB(col - 1, row));
//            rgbArrNext = rgbIntToArray(pic.getRGB(col + 1, row));
//            rgbArrPrev = rgbIntToArray(pixels[col - 1][row]);
            rgbArrPrev = rgbIntToArray(pixels[row][col - 1]);
//            rgbArrNext = rgbIntToArray(pixels[col + 1][row]);
            rgbArrNext = rgbIntToArray(pixels[row][col + 1]);
        }
        else{
//            rgbArrPrev = rgbIntToArray(pic.getRGB(col, row-1));
//            rgbArrNext = rgbIntToArray(pic.getRGB(col, row+1));
//            rgbArrPrev = rgbIntToArray(pixels[col][row-1]);
            rgbArrPrev = rgbIntToArray(pixels[row-1][col]);
//            rgbArrNext = rgbIntToArray(pixels[col][row+1]);
            rgbArrNext = rgbIntToArray(pixels[row+1][col]);
        }

        assert(rgbArrCurr.length == 3);
        for(int i = 0; i < 3; i++){
            ex += Math.pow((rgbArrNext[i] - rgbArrPrev[i]),2);
        }
        return ex;
    }

    private Queue<Integer> adjVerticalVertex(int v){
        Queue<Integer> q = new Queue<>();
        // If we are dealing with the virtual source, all pixels on the first row are adjacent to it.
        if(v==srcVertex){
            int row = 0;
            for(int col = 0; col < W; col++) {
                assert isValidPixel(col,row);
                int w = colrowToVertex(col,row);
                q.enqueue(w);
            }
            return q;
        }

        // If we are dealing with sink, return empty q.
        if(v == sinkVertex)
            return q;

        int[] colRow = vertexToColRow(v);
        int currCol = colRow[0];
        int currRow = colRow[1];
        assert isValidPixel(currCol,currRow);

        // If we are dealing with last row, then adjacent vertex is sink source.
        if(currRow == H-1){
            q.enqueue(sinkVertex);
            return q;
        }

        int[] neighborColOffset = {-1,0,+1};
        for(int offset : neighborColOffset){
            int col = currCol + offset;
            int row = currRow + 1;
            if(isValidPixel(col,row))
                q.enqueue(colrowToVertex(col,row));
        }
        return q;
    }

    private double weight(int v,int w){
        // edge weight is defined by energy at node w.
        // for w = sink vertex, define weight to be 0.
        assert(w != srcVertex);
        if(w == sinkVertex)
            return 0.0;

        int[] colRow = vertexToColRow(w);
        return energy(colRow[0],colRow[1]);
    }

    private boolean isValidPixel(int col, int row){
        if(col < 0 || col >= W || row < 0 || row >= H)
            return false;
        else
            return true;
    }

    private boolean isValidVertex(int v){
        return v < nVertices();
    }

    private void transpose(){
        int H = this.H;
        this.H = this.W;
        this.W = H;
        pixels = transpose(pixels);
        pixEnergyArr = transpose(pixEnergyArr);
        isTransposed = !isTransposed;
    }
    private int[][] transpose(int[][] A){
        int m, n;
//        if(!isTransposed){
            m = this.H;
            n = this.W;
//        }
//        else{
//            m = this.W;
//            n = this.H;
//        }
        int[][] B = new int[n][m];

        for(int row = 0; row < m; row++){
            for(int col = 0; col < n; col++){
                B[col][row] = A[row][col];
            }
        }
        return B;
    }
    private double[][] transpose(double[][] A){
        int m, n;
//        if(!isTransposed){
            m = this.H;
            n = this.W;
//        }
//        else{
//            m = this.W;
//            n = this.H;
//        }
        double[][] B = new double[n][m];

        for(int row = 0; row < m; row++){
            for(int col = 0; col < n; col++){
                B[col][row] = A[row][col];
            }
        }
        return B;
    }
    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture){
        pic = new Picture(picture);

        this.W = pic.width();
        this.H = pic.height();
        int width = this.W;
        int height = this.H;

        // Store local copy of pixel values
//        pixels = new int[W][H];
        pixels = new int[H][W];
        pixEnergyArr2 = new double[H][W];
        for(int col = 0; col < W; col++){
            for(int row = 0; row < H; row++){
//                pixels[col][row] = pic.getRGB(col,row);
                pixels[row][col] = pic.getRGB(col,row);
            }
        }
        for(int col = 0; col < W; col++){
            for(int row = 0; row < H; row++){
                pixEnergyArr2[row][col] = computePixEnergy(col,row);
            }
        }

//        pixEnergy = new double[width][height];
        pixEnergyArr = new double[height][width];
        // Store border pixel energy
        // Top and bottom rows
        int[] borderRows = {0, height -1};
        for(int row : borderRows){
            for(int col = 0; col < width; col++){
//                pixEnergy[col][row] = BORDER_ENERGY;
                pixEnergyArr[row][col] = BORDER_ENERGY;
            }
        }
        // Left and right columns.
        int[] borderCols = {0, width -1};
        for(int col : borderCols){
            for(int row = 0; row < height; row++){
//                pixEnergy[col][row] = BORDER_ENERGY;
                pixEnergyArr[row][col] = BORDER_ENERGY;
            }
        }

        // Compute interior pixel energy
        for(int col = 1; col < width -1; col++){
            for(int row = 1; row < height -1; row++){
                double ex = energyCol(col, row,false);
                double ey = energyCol(col, row,true);
//                pixEnergy[col][row] = Math.sqrt(ex+ey);
                pixEnergyArr[row][col]= Math.sqrt(ex+ey);
            }
        }

        boolean energyCheck = Arrays.deepEquals(pixEnergyArr2,pixEnergyArr);
        assert energyCheck;


        // Define the weight of an edge v->w as the energy of pixel w
        // We create a virtual vertex above and below picture as sources and sink.

        int nVertex = nPixels()+2;
        edgeTo = new int[nVertex];
        distTo = new double[nVertex];
        for(int v = 0; v < nVertex; v++)
            distTo[v] = Double.POSITIVE_INFINITY;
        distTo[0] = 0.0;
    }

    // current picture
    public Picture picture(){
        if(seamRemoved){
            pic = new Picture(this.W, this.H);
            for (int col = 0; col < this.W; col++) {
                for (int row = 0; row < this.H; row++) {
                    if(isTransposed)
                        pic.setRGB(col,row,pixels[row][col]);
                    else
                        pic.setRGB(col,row,pixels[col][row]);
                }
            }
            seamRemoved = false;
        }
        return pic;
    }

    // width of current picture
    public int width(){
        if(isTransposed)
            return this.H;
        else
            return this.W;
    }

    // height of current picture
    public int height(){
        if(isTransposed)
            return this.W;
        else
            return this.H;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y){
        return pixEnergyArr[y][x];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam(){
        if(!isTransposed){
            transpose();
//            pixels = transpose(pixels);
//            pixEnergy = transpose(pixEnergy);
//            isTransposed = true;
        }
        return findVerticalSeam();
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam(){
        if(isTransposed){
            transpose();
        }
        // Topological order = put vertices in order such that all its directed edges point from a vertex earlier in order to vertex later in order. For a vertical seam in a picture, that order is by traversing pixels row by row
        int[] edgeTo = new int[nVertices()];
        double[] distTo = new double[nVertices()];
        Arrays.fill(distTo,Double.POSITIVE_INFINITY);
        distTo[srcVertex] = 0.0;
        for(int v = 0; v < nVertices(); v++){
            for(int w : adjVerticalVertex(v)){
                double currWeight = weight(v,w);
                if(distTo[w] >= distTo[v] + currWeight){
                    distTo[w] = distTo[v] + currWeight;
                    edgeTo[w] = v;
                }
            }
        }

        Stack<Integer> s = new Stack<>();
        for(int v = edgeTo[sinkVertex]; v != srcVertex; v = edgeTo[v]){
            s.push(v);
        }

        int[] seam = new int[s.size()];
        for(int i = 0; i < s.size(); i++) {
            int[] colRow = vertexToColRow(s.pop());
            // Convert vertex to corresponding col,row and store column
            seam[i] = colRow[0];
        }

        return seam;
    }

    // remove horizontal seam from current picture
//    public void removeHorizontalSeam(int[] seam){}

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam){
        if(isTransposed){
            transpose();
//            pixels = transpose(pixels);
//            pixEnergy = transpose(pixEnergy);
//            isTransposed = false;
        }
//        int srcPos;
//        int destPos;
//        int length;
        int[][] pixelsNew = new int[H][W-1];
        double[][] energyNew = new double[H][W-1];
        int rowCnt = 0;
        for(int col : seam){
           if(col == 0){
               int srcPos = 1;
               int destPos = 0;
               int length = W-1;
               System.arraycopy(pixels[rowCnt],srcPos,pixelsNew[rowCnt],destPos,length);
               System.arraycopy(pixEnergyArr[rowCnt],srcPos,energyNew[rowCnt],destPos,length);
           }
           else if (col == W-1) {
               int srcPos = 0;
               int destPos = 0;
               int length = W-1;
               System.arraycopy(pixels[rowCnt],srcPos,pixelsNew[rowCnt],destPos,length);
               System.arraycopy(pixEnergyArr[rowCnt],srcPos,energyNew[rowCnt],destPos,length);
           }
           else{
               int[] srcPosAll = {0, col+1};
               int[] destPosAll = {0, col};
               int[] lengthAll = {col, W-1-col};
               for(int i = 0; i < srcPosAll.length; i++){
                   int srcPos = srcPosAll[i];
                   int destPos =  destPosAll[i];
                   int length = lengthAll[i];
                   System.arraycopy(pixels[rowCnt],srcPos,pixelsNew[rowCnt],destPos,length);
                   System.arraycopy(pixEnergyArr[rowCnt],srcPos,energyNew[rowCnt],destPos,length);
               }
           }
            rowCnt++;
        }
        pixels = pixelsNew;
        pixEnergyArr = energyNew;
        seamRemoved = true;
        W--;
    }

    // unit testing
    public static void main(String[] args){
        Picture picture = new Picture(args[0]);
        SeamCarver sc = new SeamCarver(picture);

        for(int row = 0; row < sc.height(); row++){
            for(int col = 0; col < sc.width(); col++){
                StdOut.printf("%9.0f",sc.energy(col,row));
            }
            StdOut.println();
        }
        StdOut.println();

        for(int row = 0; row < sc.height(); row++){
            for(int col = 0; col < sc.width(); col++){
                StdOut.printf("%9.0f",sc.pixEnergyArr2[row][col]);
            }
            StdOut.println();
        }
    }
}