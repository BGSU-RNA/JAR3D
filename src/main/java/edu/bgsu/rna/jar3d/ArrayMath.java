package edu.bgsu.rna.jar3d;
public class ArrayMath {
	public static double mean(double[] arr){
		int n = arr.length;
		double sum = 0;
		for(int i = 0; i < n; i++)
		{
			sum = sum + arr[i];
		}
		double mean = (double)sum / (double)n;
		return mean;
	}
	public static double min(double[] arr){
		int n = arr.length;
		double min = arr[0];
		for(int i = 0; i < n; i++)
		{
			min = Math.min(min, arr[i]);
		}
		return min;
	}
	public static double max(double[] arr){
		int n = arr.length;
		double max = arr[0];
		for(int i = 0; i < n; i++)
		{
			max = Math.max(max, arr[i]);
		}
		return max;
	}
	public static double median(double[] arr){
		if(arr.length == 1) {
			return (arr[0]);
		}
		int n = arr.length;
		java.util.Arrays.sort(arr);
		double middle = ((double)n)/2;  // subscript of middle element
	    if (n%2 == 1) {
	        // Odd number of elements -- return the middle one.
	        return arr[(int)middle];
	    } else {
	       // Even number -- return average of middle two
	       // Must cast the numbers to double before dividing.
	       return (arr[(int)middle] + arr[(int)(middle)-1]) / 2.0;
	    }
	}
	public static double mean(int[] arr){
		int n = arr.length;
		double sum = 0;
		for(int i = 0; i < n; i++)
		{
			sum = sum + arr[i];
		}
		double mean = (double)sum / (double)n;
		return mean;
	}
	public static int min(int[] arr){
		int n = arr.length;
		int min = arr[0];
		for(int i = 0; i < n; i++)
		{
			min = Math.min(min, arr[i]);
		}
		return min;
	}
	public static int max(int[] arr){
		int n = arr.length;
		int max = arr[0];
		for(int i = 0; i < n; i++)
		{
			max = Math.max(max, arr[i]);
		}
		return max;
	}
	public static double median(int[] arr){
		if(arr.length == 1) return (double)(arr[0]);
		int n = arr.length;
		java.util.Arrays.sort(arr);
		double middle = ((double)n)/2;  // subscript of middle element
	    if (n%2 == 1) {
	        // Odd number of elements -- return the middle one.
	        return arr[(int)middle];
	    } else {
	       // Even number -- return average of middle two
	       // Must cast the numbers to double before dividing.
	       return (arr[(int)middle] + arr[(int)(middle)-1]) / 2.0;
	    }
	}
}
