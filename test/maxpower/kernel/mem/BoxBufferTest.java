package maxpower.kernel.mem;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import com.maxeler.maxcompiler.v2.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.managers.SimulationParams;
import com.maxeler.maxcompiler.v2.managers.standard.SimulationManager;
import com.maxeler.maxcompiler.v2.utils.Bits;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

public class BoxBufferTest {

	private static final int m_itemBitWidth = 24;

	@Test public void coprime_5_3()          { testBoxBuffer(1, 3125,  5,  3); }
	@Test public void simple_2_2()           { testBoxBuffer(1, 1024,  2,  2); }
	@Test public void nonPow2simle()         { testBoxBuffer(1, 2508, 12, 12); }
	@Test public void pow2multiple_4_16()    { testBoxBuffer(1, 2500,  4, 16); }
	@Test public void nonPow2multiple_4_20() { testBoxBuffer(1, 2500,  4, 20); }
	@Test public void nonMultiple_12_16()    { testBoxBuffer(1, 2508, 12, 16); }
	@Test public void nonMultiple_16_20()    { testBoxBuffer(1, 2512, 16, 20); }
	@Test public void pow2factor_16_8()      { testBoxBuffer(1, 2512, 16,  8); }
	@Test public void nonFactor_12_8()       { testBoxBuffer(1, 2508, 12,  8); }

	@Test public void simple_2d()            { testBoxBuffer(2, 64, 2, 2, 2); }
	@Test public void coprime_2d()           { testBoxBuffer(2, 90, 5, 3, 2); }
	@Test public void coprime_3d()           { testBoxBuffer(3, 90, 5, 3, 2, 5); }
	@Test public void coprime_4d()           { testBoxBuffer(4, 24, 4, 3, 3, 3, 3); }


	private void testBoxBuffer(int numDimensions, int maxItemsPerDim, int numInputItems, int... numOutputItems) {
		SimulationManager mgr = new SimulationManager("BoxBufferTest_"+numDimensions+"D_"+maxItemsPerDim+"_"+numInputItems+"_"+numOutputItems[0],  SimulationParams.BITACCURATE_MAX4);
		int[] maxItems = new int[numDimensions];
		for (int i = 0; i < numDimensions; i++) {
			maxItems[i] = maxItemsPerDim;
		}

		TestData data = new TestData(maxItems, numInputItems, numOutputItems);
		TestKernel dutA = new TestKernel(mgr.makeKernelParameters(), numDimensions, maxItems, numInputItems, numOutputItems);

		mgr.setKernel(dutA);

		mgr.setInputDataRaw("wrData", data.m_wrData);

		mgr.setInputData("wrEnable", data.m_wrEnable);
		mgr.setInputData("wrBuffer", data.m_wrBuffer);
		mgr.setInputData("rdBuffer", data.m_rdBuffer);
		for (int i = 0; i < numDimensions; i++) {
			mgr.setInputData("wrIndex"+i,  data.m_wrIndex[i]);
			mgr.setInputData("rdIndex"+i,  data.m_rdIndex[i]);
		}

		mgr.setKernelCycles(data.m_numCycles);

		mgr.runTest();

		assertTrue(data.testOutput(mgr.getOutputDataRaw("rdData")));
	}


	private class TestKernel extends Kernel {
		TestKernel(KernelParameters p, int numDimensions, int[] maxItems, int numInputItems, int[] numOutputItems) {
			super(p);
			int numCycles = product(maxItems);
			DFEVectorType<DFEVar> inType  = new DFEVectorType<DFEVar>(dfeUInt(m_itemBitWidth), numInputItems);
			DFEVectorType<DFEVar> outType = new DFEVectorType<DFEVar>(dfeUInt(m_itemBitWidth), product(numOutputItems));
			DFEVar[] wrIndex = new DFEVar[numDimensions];
			DFEVar[] rdIndex = new DFEVar[numDimensions];
			for (int i = 0; i < numDimensions; i++) {
				int idxBits = MathUtils.bitsToAddress(maxItems[i]);
				wrIndex[i] = io.input("wrIndex"+i,  dfeUInt(idxBits)).simWatch("wrIndex"+i);
				rdIndex[i] = stream.offset(io.input("rdIndex"+i,  dfeUInt(idxBits)).simWatch("rdIndex"+i), -numCycles);
			}

			DFEVar wrBuffer = io.input("wrBuffer", dfeUInt(1));
			DFEVar wrEnable = io.input("wrEnable", dfeUInt(1));

			DFEVar rdBuffer = stream.offset(io.input("rdBuffer", dfeUInt(1)), -numCycles);

			DFEVector<DFEVar> wrData = io.input("wrData", inType);
			wrData.simWatch("wrData");

			BoxBuffer<DFEVar> buffer = new BoxBuffer<DFEVar>(this, maxItems, numOutputItems, inType);
			buffer.write(wrData, wrIndex, wrEnable, wrBuffer);
			io.output("rdData", outType) <== stream.offset(buffer.read(rdIndex, rdBuffer), numCycles).simWatch("rdData");
		}
	}

	private class TestData {
		final int[] m_data;
		final Bits[] m_wrData;
		final double[] m_wrEnable;
		final double[] m_wrBuffer;
		final double[] m_rdBuffer;
		final double[][] m_wrIndex;
		final double[][] m_rdIndex;
		final int m_numCycles;
		final int[] m_numOutputItems;
		final int[] m_maxItems;


		TestData(int[] maxItems, int numInputItems, int[] numOutputItems) {
			DFEVectorType<DFEVar> inType = new DFEVectorType<DFEVar>(Kernel.dfeUInt(m_itemBitWidth), numInputItems);
			if (maxItems[maxItems.length - 1] % numInputItems != 0) {
				throw new RuntimeException("Maximum number of items in the fast dimension needs to be a multiple of the number of input items.");
			}
			m_numCycles       = product(maxItems) / numInputItems;
			m_numOutputItems  = new int[numOutputItems.length];
			m_maxItems        = new int[maxItems.length];
			int[] lastAddress = new int[maxItems.length];
			for (int i = 0; i < maxItems.length; i++) {
				m_numOutputItems[i] = numOutputItems[i];
				m_maxItems[i] = maxItems[i];
				lastAddress[i] = maxItems[i] - numOutputItems[i];
			}
			int lastIndex = getIndex(lastAddress);

			m_data = new int[product(m_maxItems)];
			for (int i = 0; i < m_data.length; i++) {
				m_data[i] = i % (1 << m_itemBitWidth);
			}

			m_wrData   = new Bits[m_numCycles];
			m_wrEnable = new double[m_numCycles];
			m_wrBuffer = new double[m_numCycles];
			m_rdBuffer = new double[m_numCycles];
			m_wrIndex  = new double[m_maxItems.length][m_numCycles];
			m_rdIndex  = new double[m_maxItems.length][m_numCycles];

			// Build up input data
			for (int i = 0; i < m_numCycles; i++) {
				m_wrEnable[i] = 1;
				m_wrBuffer[i] = 0;
				m_rdBuffer[i] = 0;
				int[] input = new int[numInputItems];
				for (int j = 0; j < numInputItems; j++) {
					input[j] = m_data[numInputItems * i + j];
				}
				m_wrData[i] = inType.encodeConstant(input);

				int readIndex = 0;
				// Test the corners first, then random offsets
				switch (i) {
					case 0:
						readIndex = lastIndex;
						break;
					case 1:
						readIndex = 0;
						break;
					default:
						readIndex = ((int) (Math.random() * lastIndex));
				}
				int[] readIndices = getAddress(readIndex);
				for (int dim = 0; dim < readIndices.length; dim++) {//Shift it back into a valid range if we have gone off the edge
					readIndices[dim] = Math.min(readIndices[dim], m_maxItems[dim] - m_numOutputItems[dim]);
				}
				int[] writeIndices = getAddress(numInputItems * (long)i);
				for (int j = maxItems.length - 1; j >= 0; j--) {
					m_wrIndex[j][i] = writeIndices[j];
					m_rdIndex[j][i] = readIndices[j];
				}
			}
		}

		int[] getAddress(long index) {
			int denom = 1;
			int[] address = new int[m_maxItems.length];
			for (int dim = m_maxItems.length - 1; dim >= 0; dim--) {
				address[dim] = (int) ((index / denom) % m_maxItems[dim]);
				denom *= m_maxItems[dim];
			}
			return address;
		}

		int getIndex(int[] address) {
			int skip = 1;
			int index = 0;
			for (int dim = m_maxItems.length - 1; dim >= 0; dim--) {
				index += address[dim] * skip;
				skip *= m_maxItems[dim];
			}
			return index;
		}

		int[] extractBlock(int[] address) {
			int[] output = new int[product(m_numOutputItems)];
			for (int i = 0; i < output.length; i++) {
				int[] newAddress = new int[address.length];
				int denom = 1;
				for (int dim = m_numOutputItems.length - 1; dim >= 0; dim--) {
					newAddress[dim] = address[dim] + (i / denom) % m_numOutputItems[dim];
					denom *= m_numOutputItems[dim];
				}
				int index = getIndex(newAddress);
				output[i] = m_data[index];
			}
			return output;
		}

		int[] getReadAddress(int cycle) {
			int[] address = new int[m_maxItems.length];
			for (int dim = 0; dim < m_maxItems.length; dim++) {
				address[dim] = (int)m_rdIndex[dim][cycle];
			}
			return address;
		}

		String printArray(int[] array) {
			String output = " { ";
			for (int i = 0; i < array.length - 1; i++) {
				output += array[i] + ", ";
			}
			output += array[array.length - 1] + " } ";
			return output;
		}

		boolean testOutput(List<Bits> rdData) {
			DFEVectorType<DFEVar> outType = new DFEVectorType<DFEVar>(Kernel.dfeUInt(m_itemBitWidth), product(m_numOutputItems));
			boolean testPassed = true;
			for (int i = 0; i < m_numCycles; i++) {
				@SuppressWarnings("unchecked")
                List<Double> output = outType.decodeConstant(rdData[i]);
				int[] expected = extractBlock(getReadAddress(i));
				for (int j = 0; j < expected.length; j++) {
					if (output[j].intValue() != expected[j]) {
						System.out.println("[" + i + ", " + j + "] " + printArray(getReadAddress(i)) + " Value expected: " + expected[j] + " != got: " + output[j].intValue());
						testPassed = false;
					}
				}
			}
			return testPassed;
		}
	}


	public static int product(int[] x) {
		int result = x[0];
		for (int i = 1; i < x.length; i++) {
			result *= x[i];
		}
		return result;
	}
}
