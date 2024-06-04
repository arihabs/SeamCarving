import edu.princeton.cs.algs4.Picture;

public class SeamCarver {
    private Picture pic;
    private double[][] pixEnergy;
    private final double BORDER_ENERGY = 1000;

    private int[] rgbIntToArray(int rgb){
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >>  8) & 0xFF;
        int b = (rgb >>  0) & 0xFF;
        int[] arr = {r,g,b};
        return arr;
    }

    private double energyCol(int col, int row, boolean transpose){
//        assert(col > 0 && col < width()-1);
//        assert(row > 0 && col < height()-1);
        double ex = 0.0;
        int[] rgbArrCurr = rgbIntToArray(pic.getRGB(col,row));
        int[] rgbArrPrev;
        int[] rgbArrNext;
        if(!transpose) {
            rgbArrPrev = rgbIntToArray(pic.getRGB(col - 1, row));
            rgbArrNext = rgbIntToArray(pic.getRGB(col + 1, row));
        }
        else{
            rgbArrPrev = rgbIntToArray(pic.getRGB(col, row-1));
            rgbArrNext = rgbIntToArray(pic.getRGB(col, row+1));
        }

        assert(rgbArrCurr.length == 3);
        for(int i = 0; i < 3; i++){
            ex += Math.pow((rgbArrNext[i] - rgbArrPrev[i]),2);
        }
        return ex;
    }
    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture){
        pic = new Picture(picture);

        pixEnergy = new double[width()][height()];
        // Store border pixel energy
        // Top and bottom rows
        int[] borderRows = {0, height()-1};
        for(int row : borderRows){
            for(int col = 0; col < width(); col++){
                pixEnergy[col][row] = BORDER_ENERGY;
            }
        }
        // Left and right columns.
        int[] borderCols = {0, width()-1};
        for(int col : borderCols){
            for(int row = 0; row < height(); row++){
                pixEnergy[col][row] = BORDER_ENERGY;
            }
        }

        // Compute interior pixel energy
        for(int col = 1; col < width()-1; col++){
            for(int row = 1; row < height()-1; row++){
                double ex = energyCol(col, row,false);
                double ey = energyCol(col, row,true);
                pixEnergy[col][row] = Math.sqrt(ex+ey);
            }
        }
    }

    // current picture
    public Picture picture(){
        return pic;
    }

    // width of current picture
    public int width(){
        return pic.width();
    }

    // height of current picture
    public int height(){
        return pic.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y){
        return pixEnergy[x][y];
    }

    // sequence of indices for horizontal seam
//    public int[] findHorizontalSeam(){}

    // sequence of indices for vertical seam
//    public int[] findVerticalSeam(){}

    // remove horizontal seam from current picture
//    public void removeHorizontalSeam(int[] seam){}

    // remove vertical seam from current picture
//    public void removeVerticalSeam(int[] seam){}

    // unit testing
    public static void main(String[] args){}
}