import edu.princeton.cs.algs4.Picture;

public class SeamCarver {
    private Picture pic;
    private double[][] energy;
    private final double BORDER_ENERGY = 1000;

    private int[] rgbIntToArray(int rgb){
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >>  8) & 0xFF;
        int b = (rgb >>  0) & 0xFF;
        int[] arr = {r,g,b};
        return arr;
    }

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture){
        pic = new Picture(picture);

        // Compute interior pixel energy
        for(int col = 1; col < width()-1; col++){
            for(int row = 1; row < height()-1; row++){
                int rgb = pic.getRGB(col,row);
                int[] rgbArr = rgbIntToArray(rgb);
            }
        }

        // Set border energy to predefinde value
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
    public double energy(int x, int y){}

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam(){}

    // sequence of indices for vertical seam
    public int[] findVerticalSeam(){}

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam){}

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam){}

    // unit testing
    public static void main(String[] args){}
}