import java.io.*;
import java.util.Arrays;

public class ReadExample
{
	public static void main(String[] args)
	{
		String PATH = "/Users/macbook/Documents/workspace/SoundWatchPlayground/src/example.wav";
		String LG_PATH = "/Users/macbook/Documents/workspace/SoundWatchPlayground/src/example_LG.wav";
		String VIDEO2 = "/Users/macbook/Documents/workspace/SoundWatchPlayground/src/example2.wav";
		try
		{
			// Open the wav file specified as the first argument
			WavFile wavFile = WavFile.openWavFile(new File(VIDEO2));

			// Display information about the wav file
			wavFile.display();

			// Get the number of audio channels in the wav file
			int numChannels = wavFile.getNumChannels();
			System.out.println("numChannels" + numChannels);
			System.out.println("Num frames: " + wavFile.getNumFrames());
			// Create a buffer of 100 frames
			double[] buffer = new double[(int) wavFile.getNumFrames() * numChannels];
//			int[][] buffer2 = new int[(int) wavFile.getNumFrames()][numChannels];
			int[][] buffer2 = new int[numChannels][(int) wavFile.getNumFrames()];
			int[] buffer3 = new int[(int) wavFile.getNumFrames() * numChannels];



			int framesRead;
			double min = Double.MAX_VALUE;
			double max = Double.MIN_VALUE;

			do
			{
				// Read frames into buffer
				framesRead = wavFile.readFrames(buffer, buffer.length);
//				framesRead = wavFile.readFrames(buffer2, 0, buffer.length);
//				framesRead = wavFile.readFrames(buffer3, buffer.length);
				// Loop through frames and look for minimum and maximum value
				for (int s=0 ; s<framesRead * numChannels ; s++)
				{
					if (buffer[s] > max) max = buffer[s];
					if (buffer[s] < min) min = buffer[s];
				}
			}
			while (framesRead != 0);

			int LIMIT = 1000;
			int k = 0;
			for (double val : buffer) {
				System.out.print(val + " , ");
				k++;
				if (k > LIMIT) {
					break;
				}
			}


			//CODE TO PRINT OUT THE BUFFER SIZE
//			System.out.println("int 2d buffer size: " + buffer2.length + ", " + buffer2[0].length);
//			int LIMIT = 1000;
//			int k = 0;
//			for (int i = 0; i < buffer2.length; i++) {
//				for (int j = 0; j < buffer2[0].length; j++) {
//					System.out.print(buffer2[i][j] + " , ");
//					k++;
//					if (k > LIMIT) {
//						break;
//					}
//				}
//				System.out.println("");
//			}




			double[] input = new double[16000];
			for (int i = 0; i < input.length; i++) {
				input[i] = (double) buffer[i + 0];
			}
//			double[] input = buffer;


			MFCC mfcc = new MFCC();
			System.out.println("Audio dimensions: " + buffer.length);
			double[][] melSpectrogram = mfcc.melSpectrogram(input);
			System.out.println("Dimension: " + melSpectrogram.length + " , " + melSpectrogram[0].length);

////			TAKING LOG
//			1st Approach: Scale it manually by making the mean of this mel spectrogram vector to be equal to python version
			System.out.println("Log Mel Spectrogram");
			double SCALE = 1.0/7.7;
			double OFFSET = 0;
			for (int i = 0; i < melSpectrogram.length; i++) {
				for (int j = 0; j < melSpectrogram[0].length; j++) {
					melSpectrogram[i][j] = Math.log(melSpectrogram[i][j] * 1 + OFFSET) * SCALE;
				}
			}

//			float[] mfccInput = mfcc.process(input);
//			System.out.println("MFCC O")
//			System.out.println(Arrays.toString(mfccInput));
//			System.out.println(mfccInput.length);

			double sum = 0;
			double count = 0;
			min = Double.MAX_VALUE;
			max = Double.MIN_VALUE;
			System.out.println("TAKING MEL SPECTROGRAM");
			for (int i = 0; i < melSpectrogram.length; i++) {
				for (int j = 0; j < melSpectrogram[0].length; j++) {
					sum += melSpectrogram[i][j];
					count++;
					min = Math.min(min, melSpectrogram[i][j]);
					max = Math.max(max, melSpectrogram[i][j]);
				}


				for (int j = 0; j < melSpectrogram[0].length; j++) {
					System.out.print(melSpectrogram[i][j] * 1 + " , ");
				}
				System.out.println("");
				if (i > 30) {
					break;
				}
			}

			System.out.println("Average: " + sum / count);
			// Close the wavFile
			wavFile.close();

			// Output the minimum and maximum value
			System.out.printf("Min: %f, Max: %f\n", min, max);
		}
		catch (Exception e)
		{
			System.err.println(e);
		}
	}
}
