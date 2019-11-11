import java.io.*;

public class ReadExample
{
	public static void main(String[] args)
	{
		String PATH = "/Users/macbook/Documents/workspace/SoundWatchPlayground/src/example.wav";
		String LG_PATH = "/Users/macbook/Documents/workspace/SoundWatchPlayground/src/example_LG.wav";
		try
		{
			// Open the wav file specified as the first argument
			WavFile wavFile = WavFile.openWavFile(new File(PATH));

			// Display information about the wav file
			wavFile.display();

			// Get the number of audio channels in the wav file
			int numChannels = wavFile.getNumChannels();
			System.out.println("numChannels" + numChannels);
			System.out.println("Num frames: " + wavFile.getNumFrames());
			// Create a buffer of 100 frames
			double[] buffer = new double[(int) wavFile.getNumFrames() * numChannels];

			int framesRead;
			double min = Double.MAX_VALUE;
			double max = Double.MIN_VALUE;

			do
			{
				// Read frames into buffer
				framesRead = wavFile.readFrames(buffer, buffer.length);

				// Loop through frames and look for minimum and maximum value
				for (int s=0 ; s<framesRead * numChannels ; s++)
				{
					if (buffer[s] > max) max = buffer[s];
					if (buffer[s] < min) min = buffer[s];
				}
			}
			while (framesRead != 0);

			MFCC mfcc = new MFCC();
			System.out.println("Audio dimensions: " + buffer.length);
			double[][] melSpectrogram = mfcc.melSpectrogram(buffer);
			System.out.println("Dimension: " + melSpectrogram.length + " , " + melSpectrogram[0].length);

			for (int i = 0; i < melSpectrogram.length; i++) {
				for (int j = 0; j < melSpectrogram[0].length; j++) {
					melSpectrogram[i][j] = Math.log(melSpectrogram[i][j]);
				}
			}

			System.out.println("TAKING MEL SPECTROGRAM");
//			for (int i = 0; i < melSpectrogram.length; i++) {
//				for (int j = 0; j < melSpectrogram[0].length; j++) {
//					System.out.print(melSpectrogram[i][j] + " , ");
//				}
//				System.out.println("");
//			}

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
