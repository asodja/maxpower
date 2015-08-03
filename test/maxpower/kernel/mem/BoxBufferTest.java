package maxpower.kernel.mem;

import static org.junit.Assert.assertTrue;

import java.util.List;

import maxpower.kernel.mem.BoxBuffer.BoxBufferParams;

import org.junit.Test;

import com.maxeler.maxcompiler.v2.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEComplex;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEComplexType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.managers.SimulationParams;
import com.maxeler.maxcompiler.v2.managers.standard.SimulationManager;
import com.maxeler.maxcompiler.v2.utils.Bits;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

public class BoxBufferTest {

	@Test public void edgeCase_5_3() { testBoxBuffer( 500,  24,  5,  3); }
	@Test public void pow2mode()     { testBoxBuffer(1024, 128,  2,  2); }
	@Test public void crs8pipe()     { testBoxBuffer(2508,  48, 12,  8); }
	@Test public void crs12pipe()    { testBoxBuffer(2508,  48, 12, 12); }
	@Test public void crs16pipe()    { testBoxBuffer(2508,  48, 12, 16); }
	@Test public void gsmp_4_16()    { testBoxBuffer(2500,  24,  4, 16); }
	@Test public void gsmp_4_20()    { testBoxBuffer(2000,  24,  4, 20); }
	@Test public void gsmp_16_20()   { testBoxBuffer(2512,  24, 16, 20); }
	@Test public void gsmp_12_8()    { testBoxBuffer(2508,  24, 12,  8); }
	@Test public void gsmp_16_8()    { testBoxBuffer(2512,  24, 16,  8); }


	private class TestKernel extends Kernel {
		private TestKernel(KernelParameters p, int maxItems, int numInputItems, int numOutputItems, int totalBits) {
			super(p);
			BoxBufferParams params = new BoxBufferParams(maxItems, numInputItems, numOutputItems, totalBits, this);

			if (totalBits % 2 != 0)
				throw new RuntimeException("Total bits must be mod 2, because we convert to DFEComplex of n/2 bits");
			DFEComplexType complexType = new DFEComplexType(dfeUInt(totalBits / 2));
			DFEVectorType<DFEComplex> wrType = new DFEVectorType<DFEComplex>(complexType, params.numInputItems);

			DFEVar wrBuffer = io.input("wrBuffer", dfeUInt(1));
			DFEVar rdBuffer = io.input("rdBuffer", dfeUInt(1));
			DFEVar wrRow    = io.input("wrRow", dfeUInt(MathUtils.bitsToAddress(params.maxItems / params.numInputItems)));
			DFEVar rdIndex  = io.input("rdIndex", dfeUInt(MathUtils.bitsToAddress(params.maxItems)));
			DFEVar wrEnable = io.input("wrEnable", dfeUInt(1));

			DFEVector<DFEComplex> wrData = io.input("wrData", wrType);

			BoxBuffer<DFEComplex> tdRam = new BoxBuffer<DFEComplex>(this);

			DFEVectorType<DFEComplex> outType = new DFEVectorType<DFEComplex>(complexType, params.numOutputItems);
			DFEVector<DFEComplex> out = tdRam.readBox(wrBuffer, wrRow, wrData, rdBuffer, rdIndex, params, wrEnable);
			io.output("rdData", outType).connect(out);
		}
	}

	private class TestData {
		final Bits[] m_wrData;
		final Bits[] m_wrRow;
		final Bits[] m_wrEnable;
		final Bits[] m_wrBuffer;
		final Bits[] m_rdBuffer;
		final double[] m_rdIndex;
		final int numCycles;
		final int m_maxItems;
		final int m_itemBitWidth;
		final int m_numInputItems;
		final int m_numOutputItems;


		private TestData(int maxItems, int itemBitWidth, int numInputItems, int numOutputItems) {
			int inDataWidth  = numInputItems * itemBitWidth;
			int rowAddrWidth = MathUtils.bitsToAddress(maxItems / numInputItems);
			numCycles        = (int) Math.ceil(2.0 * maxItems / numInputItems);
			m_maxItems       = maxItems;
			m_itemBitWidth   = itemBitWidth;
			m_numInputItems  = numInputItems;
			m_numOutputItems = numOutputItems;

			m_wrData   = new Bits[numCycles];
			m_wrRow    = new Bits[numCycles];
			m_wrEnable = new Bits[numCycles];
			m_wrBuffer = new Bits[numCycles];
			m_rdBuffer = new Bits[numCycles];
			m_rdIndex  = new double[numCycles];
			// Build up input data
			for (int i = 0; i < 2 * maxItems; i += numInputItems) {
				int x = i / numInputItems;
				m_wrData[x]   = new Bits(inDataWidth);

				m_wrBuffer[x] = new Bits(1);
				m_rdBuffer[x] = new Bits(1);

				// Write first [maxItems] items into the buffer
				// Don't over-write as it will wrap
				if ((i + numInputItems) <= maxItems) {
					m_wrEnable[x] = new Bits(1);
					m_wrEnable[x].setBits(1);

					for (int j = 0; j < numInputItems; j++) {
						m_wrData[x].setBits(j * itemBitWidth, itemBitWidth, i + j);
					}

					// Switch buffers as we reach maxItems written to buffer
					m_wrBuffer[x].setBits(0);
					m_rdBuffer[x].setBits(1);
				} else {

					m_wrEnable[x] = Bits.allZeros(1);

					for (int j = 0; j < numInputItems; j++) {
						m_wrData[x].setOthers(0);
					}

					// Switch buffers as we reach maxItems written to buffer
					m_wrBuffer[x].setBits(1);
					m_rdBuffer[x].setBits(0);
				}

				m_wrRow[x] = new Bits(rowAddrWidth);
				m_wrRow[x].setBits(x % (maxItems / numInputItems));

				// Test the corners first, then random offsets
				switch (x - maxItems / numInputItems) {
					case 0:
						m_rdIndex[x] = maxItems - numOutputItems;
						break;
					case 1:
						m_rdIndex[x] = 0;
						break;
					default:
						m_rdIndex[x] = ((long) (Math.random() * (numInputItems * (maxItems / numInputItems) - numOutputItems)));
				}
			}
		}

		boolean testOutput(List<Bits> rdData) {
			boolean testPassed = true;
			for (int i = m_maxItems / m_numInputItems; i < (2 * m_maxItems / m_numInputItems); i++) {
				Bits tmp = new Bits(m_numOutputItems * m_itemBitWidth);
				long index = (long) m_rdIndex[i];
				for (int j = 0; j < m_numOutputItems; j++) {
					tmp.setBits(j * m_itemBitWidth, m_itemBitWidth, index++);
				}

				String tmpS = tmp.valueAsHexString();
				String rdS = rdData[i].valueAsHexString();
				if (!tmp.equals(rdData[i])) {
					System.out.println("[" + i + "] " + m_rdIndex[i] + " Value expected: " + tmpS + " != got: " + rdS);
					testPassed = false;
				}
			}
			return testPassed;
		}
	}


	private void testBoxBuffer(int maxItems, int itemBitWidth, int numInputItems, int numOutputItems) {
		SimulationManager mgr = new SimulationManager("BoxBufferTest_"+maxItems+"_"+itemBitWidth+"_"+numInputItems+"_"+numOutputItems,
						                              SimulationParams.BITACCURATE_MAX4);

		TestKernel dutA = new TestKernel(mgr.makeKernelParameters(), maxItems, numInputItems, numOutputItems, itemBitWidth);

		mgr.setKernel(dutA);

		TestData data = new TestData(maxItems, itemBitWidth, numInputItems, numOutputItems);
		mgr.setInputDataRaw("wrData",   data.m_wrData);
		mgr.setInputDataRaw("wrRow",    data.m_wrRow);
		mgr.setInputDataRaw("wrEnable", data.m_wrEnable);
		mgr.setInputDataRaw("wrBuffer", data.m_wrBuffer);
		mgr.setInputDataRaw("rdBuffer", data.m_rdBuffer);
		mgr.setInputData(   "rdIndex",  data.m_rdIndex);

		mgr.setKernelCycles(data.numCycles);

		mgr.runTest();

		assertTrue(data.testOutput(mgr.getOutputDataRaw("rdData")));
	}

}
