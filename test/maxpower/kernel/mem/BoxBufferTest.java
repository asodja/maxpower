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


	private class TestKernel extends Kernel {
		private TestKernel(KernelParameters p, int maxItems, int numInputItems, int numOutputItems, int totalBits) {
			super(p);
			DFEVectorType<DFEVar> inType  = new DFEVectorType<DFEVar>(dfeUInt(totalBits), numInputItems);
			DFEVectorType<DFEVar> outType = new DFEVectorType<DFEVar>(dfeUInt(totalBits), numOutputItems);

			int rowBits = MathUtils.bitsToAddress(maxItems / numInputItems);
			int idxBits = MathUtils.bitsToAddress(maxItems);

			DFEVar wrBuffer = io.input("wrBuffer", dfeUInt(1));
			DFEVar rdBuffer = io.input("rdBuffer", dfeUInt(1));
			DFEVar wrRow    = io.input("wrRow",    dfeUInt(rowBits));
			DFEVar rdIndex  = io.input("rdIndex",  dfeUInt(idxBits));
			DFEVar wrEnable = io.input("wrEnable", dfeUInt(1));

//int numBits = 29;
//DFEVar count = control.count.makeCounterChain().addCounter(1536 / numBits, 0);
//DFEVar in = io.input("myInput", dfeRawBits(1536), count === 0);
//DFEVar state = dfeRawBits(1536).newInstance(this);
//state <== count === 0 ? in : stream.offset(state, -1) >> numBits;
//DFEVar output = state.slice(0, numBits);

			DFEVector<DFEVar> wrData = io.input("wrData", inType);

			BoxBuffer<DFEVar> buffer = new BoxBuffer<DFEVar>(this, maxItems, numOutputItems, inType);
			buffer.write(wrData, wrRow, wrEnable, wrBuffer);
			io.output("rdData", outType) <== buffer.read(rdIndex, rdBuffer);
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
		final int m_itemBitWidth;
		final int m_numOutputItems;


		private TestData(int maxItems, int itemBitWidth, int numInputItems, int numOutputItems) {
			int inDataWidth  = numInputItems * itemBitWidth;
			int rowAddrWidth = MathUtils.bitsToAddress(maxItems / numInputItems);
			numCycles        = 2 * MathUtils.ceilDivide(maxItems, numInputItems);
			m_itemBitWidth   = itemBitWidth;
			m_numOutputItems = numOutputItems;

			m_wrData   = new Bits[numCycles];
			m_wrRow    = new Bits[numCycles];
			m_wrEnable = new Bits[numCycles];
			m_wrBuffer = new Bits[numCycles];
			m_rdBuffer = new Bits[numCycles];
			m_rdIndex  = new double[numCycles];

			// Build up input data
			for (int i = 0; i < numCycles; i++) {
				// Write first [maxItems] items into the buffer
				// Don't over-write as it will wrap
				if ((i * numInputItems + numInputItems) <= maxItems) {
					m_wrEnable[i] = Bits.allOnes(1);

					m_wrData[i] = new Bits(inDataWidth);
					for (int j = 0; j < numInputItems; j++) {
						m_wrData[i].setBits(j * itemBitWidth, itemBitWidth, i * numInputItems + j);
					}

					// Switch buffers as we reach maxItems written to buffer
					m_wrBuffer[i] = Bits.allZeros(1);
					m_rdBuffer[i] = Bits.allOnes(1);
				} else {
					m_wrEnable[i] = Bits.allZeros(1);

					m_wrData[i] = Bits.allZeros(inDataWidth);

					// Switch buffers as we reach maxItems written to buffer
					m_wrBuffer[i] = Bits.allOnes(1);
					m_rdBuffer[i] = Bits.allZeros(1);
				}

				m_wrRow[i] = new Bits(rowAddrWidth);
				m_wrRow[i].setBits(i % (maxItems / numInputItems));

				// Test the corners first, then random offsets
				switch (i - numCycles / 2) {
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
			boolean testPassed = true;
			for (int i = numCycles / 2; i < numCycles; i++) {
				//TODO: generate golden output properly (actually base it on the input), and then compare integers (not Bits).
				Bits tmp = new Bits(m_numOutputItems * m_itemBitWidth);
				long index = (long) m_rdIndex[i];
				for (int j = 0; j < m_numOutputItems; j++) {
					tmp.setBits(j * m_itemBitWidth, m_itemBitWidth, index++);
				}

				if (!tmp.equals(rdData[i])) {
					System.out.println("[" + i + "] " + m_rdIndex[i] + " Value expected: " + tmp.valueAsHexString() + " != got: " + rdData[i].valueAsHexString());
					testPassed = false;
				}
			}
			return testPassed;
		}
	}
}
