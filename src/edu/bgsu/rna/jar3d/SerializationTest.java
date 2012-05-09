package edu.bgsu.rna.jar3d;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Random;

public class SerializationTest {
	
	public static void writeobj(String filename){
		Random generator = new Random();
		double[][] object = new double[10][10];
		for(int i = 0; i<10; i++){
			for(int j = 0; j<10; j++){
				object[i][j] = generator.nextDouble();
			}
		}
		FileOutputStream fos = null;
		ObjectOutputStream out = null;
		try {
		    fos = new FileOutputStream(filename);
		    out = new ObjectOutputStream(fos);
		    out.writeObject(object);
		    out.close();
		}
		catch(IOException ex) {
			ex.printStackTrace();
		}
	}
	public static void readobj(String filename){
		double[][] object = null;
		FileInputStream fis = null;
		ObjectInputStream in = null;
		try
		{
		  fis = new FileInputStream(filename);
		  in = new ObjectInputStream(fis);
		  object = (double[][])in.readObject();
		  in.close();
		}
		catch(IOException ex)
		{
		  ex.printStackTrace();
		}
		catch(ClassNotFoundException ex)
		{
		  ex.printStackTrace();
		}
		for(int i = 0; i<10; i++){
			for(int j = 0; j<10; j++){
				System.out.print(object[i][j]);
				System.out.print(" ");
			}
			System.out.println("");
		}
	}
}