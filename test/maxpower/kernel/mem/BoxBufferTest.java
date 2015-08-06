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

	@Test public void edgeCase_5_3() { testBoxBuffer( 500, 24,  5,  3); }
	@Test public void pow2mode()     { testBoxBuffer(1024, 64,  2,  2); }
	@Test public void crs8pipe()     { testBoxBuffer(2508, 48, 12,  8); }
	@Test public void crs12pipe()    { testBoxBuffer(2508, 48, 12, 12); }
	@Test public void crs16pipe()    { testBoxBuffer(2508, 48, 12, 16); }
	@Test public void gsmp_4_16()    { testBoxBuffer(2500, 24,  4, 16); }
	@Test public void gsmp_4_20()    { testBoxBuffer(2000, 24,  4, 20); }
	@Test public void gsmp_16_20()   { testBoxBuffer(2512, 24, 16, 20); }
	@Test public void gsmp_12_8()    { testBoxBuffer(2508, 24, 12,  8); }
	@Test public void gsmp_16_8()    { testBoxBuffer(2512, 24, 16,  8); }


	private void testBoxBuffer(int maxItems, int itemBitWidth, int numInputItems, int numOutputItems) {
		SimulationManager mgr = new SimulationManager("BoxBufferTest_"+maxItems+"_"+itemBitWidth+"_"+numInputItems+"_"+numOutputItems,
						                              SimulationParams.BITACCURATE_MAX4);

		TestKernel dutA = new TestKernel(mgr.makeKernelParameters(), maxItems, numInputItems, numOutputItems, itemBitWidth);

		mgr.setKernel(dutA);

		TestData data = new TestData(maxItems, itemBitWidth, numInputItems, numOutputItems);
		mgr.setInputDataRaw("wrData", data.m_wrData);

		mgr.setInputData("wrEnable", data.m_wrEnable);
		mgr.setInputData("wrBuffer", data.m_wrBuffer);
		mgr.setInputData("rdBuffer", data.m_rdBuffer);
		mgr.setInputData("wrIndex",  data.m_wrIndex);
		mgr.setInputData("rdIndex",  data.m_rdIndex);

		mgr.setKernelCycles(data.m_numCycles);

		mgr.runTest();

		assertTrue(data.testOutput(mgr.getOutputDataRaw("rdData")));
	}


	private class TestKernel extends Kernel {
		private TestKernel(KernelParameters p, int maxItems, int numInputItems, int numOutputItems, int totalBits) {
			super(p);
			DFEVectorType<DFEVar> inType  = new DFEVectorType<DFEVar>(dfeUInt(totalBits), numInputItems);
			DFEVectorType<DFEVar> outType = new DFEVectorType<DFEVar>(dfeUInt(totalBits), numOutputItems);

			int idxBits = MathUtils.bitsToAddress(maxItems);

			DFEVar wrBuffer = io.input("wrBuffer", dfeUInt(1));
			DFEVar rdBuffer = io.input("rdBuffer", dfeUInt(1));
			DFEVar wrIndex  = io.input("wrIndex",  dfeUInt(idxBits));
			DFEVar rdIndex  = io.input("rdIndex",  dfeUInt(idxBits));
			DFEVar wrEnable = io.input("wrEnable", dfeUInt(1));

			DFEVector<DFEVar> wrData = io.input("wrData", inType);

			BoxBuffer<DFEVar> buffer = new BoxBuffer<DFEVar>(this, maxItems, numOutputItems, inType);
			buffer.write(wrData, wrIndex, wrEnable, wrBuffer);
			io.output("rdData", outType) <== buffer.read(rdIndex, rdBuffer);
		}
	}

	private class TestData {
		final int[] m_data;
		final Bits[] m_wrData;
		final double[] m_wrEnable;
		final double[] m_wrBuffer;
		final double[] m_rdBuffer;
		final double[] m_wrIndex;
		final double[] m_rdIndex;
		final int m_numCycles;
		final int m_itemBitWidth;
		final int m_numOutputItems;


		private TestData(int maxItems, int itemBitWidth, int numInputItems, int numOutputItems) {//TODO: multidim
			DFEVectorType<DFEVar> inType = new DFEVectorType<DFEVar>(Kernel.dfeUInt(itemBitWidth), numInputItems);
			m_numCycles        = 2 * MathUtils.ceilDivide(maxItems, numInputItems);
			m_itemBitWidth   = itemBitWidth;
			m_numOutputItems = numOutputItems;

			m_data = new int[maxItems];
			int maxValue = itemBitWidth < 31 ? 1 << itemBitWidth : Integer.MAX_VALUE;
			for (int i = 0; i < maxItems; i++) {
				m_data[i] = i % maxValue;
			}

			m_wrData   = new Bits[m_numCycles];
			m_wrEnable = new double[m_numCycles];
			m_wrBuffer = new double[m_numCycles];
			m_rdBuffer = new double[m_numCycles];
			m_wrIndex  = new double[m_numCycles];
			m_rdIndex  = new double[m_numCycles];

			// Build up input data
			for (int i = 0; i < m_numCycles; i++) {
				boolean phase1 = i < m_numCycles / 2;
				m_wrEnable[i] = phase1 ? 1 : 0;
				m_wrBuffer[i] = phase1 ? 0 : 1;
				m_rdBuffer[i] = phase1 ? 1 : 0;
				int[] input = new int[numInputItems];
				for (int j = 0; j < numInputItems; j++) {
					input[j] = phase1 ? m_data[numInputItems * i + j] : 0;
				}
				m_wrData[i] = inType.encodeConstant(input);

				m_wrIndex[i] = (numInputItems * i) % maxItems;

				// Test the corners first, then random offsets
				switch (i - m_numCycles / 2) {
					case 0:
						m_rdIndex[i] = maxItems - numOutputItems;
						break;
					case 1:
						m_rdIndex[i] = 0;
						break;
					default:
						m_rdIndex[i] = ((long) (Math.random() * (numInputItems * (maxItems / numInputItems) - numOutputItems)));
				}
			}
		}

		boolean testOutput(List<Bits> rdData) {
			DFEVectorType<DFEVar> outType = new DFEVectorType<DFEVar>(Kernel.dfeUInt(m_itemBitWidth), m_numOutputItems);
			boolean testPassed = true;
			for (int i = m_numCycles / 2; i < m_numCycles; i++) {
				@SuppressWarnings("unchecked")
                List<Double> output = outType.decodeConstant(rdData[i]);
				for (int j = 0; j < m_numOutputItems; j++) {
					if (output[j].intValue() != m_data[(int)m_rdIndex[i] + j]) {
						System.out.println("[" + i + "] " + ((int)m_rdIndex[i] + j) + " Value expected: " + m_data[(int)m_rdIndex[i] + j] + " != got: " + output[j].intValue());
						testPassed = false;
					}
				}
			}
			return testPassed;
		}
	}
}
