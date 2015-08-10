package maxpower.kernel.mem;

import java.util.ArrayList;
import java.util.List;

import maxpower.kernel.arithmetic.ConstDenominator;
import maxpower.kernel.arithmetic.ConstDenominator.ConstDivModResult;

import com.maxeler.maxcompiler.v2.errors.MaxCompilerAPIError;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelLib;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.memory.Memory;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.KernelObjectVectorizable;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.KernelType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

/**
 * The Box buffer reads in a certain sized stream of data and holds it in a
 * buffer, from which you can query it to retrieve a box of adjacent items at
 * a random location. In this case box is meant in the geometric sense (also
 * called an n-orthotope or hyperectangle), so a contiguous section of an
 * N-dimensional space can be retrieved. In 1D this is a simple interval.
 *
 * Internally it can be double buffered, so we can write N items a cycle into the
 * buffer and read M items out. Obviously for maximum performance, we need to set-up
 * our kernel so that N*cycles >= items_to_write otherwise we can't read and write
 * at full speed. If it is not double buffered, then it is up to the user to make
 * sure that data has not been overwritten before it is read.
 *
 * If you have a section of data that you need to read chunks from which will then be
 * swapped out for another section then you should use the double buffered version so
 * that the new section can be written while the old section is read, so we can do
 * useful compute on every cycle. If instead you need a sliding window or you only
 * have a single section to read, then use single buffering.
 *
 * Note that data must be written in linearly with new data on every cycle (when buffer
 * is in use), but can be read out in any order. This is due to a stream offset used
 * to cluster input items together to use FMem more efficiently. To disable this
 * behaviour, set defaultLutCost to zero (this forces it to only optimise for LUTs).
 *
 * Also note that when using this buffer with more than one dimension, the fast dimension
 * (the one we write in linearly) is the last dimension specified.
 */
public class BoxBuffer<T extends KernelObjectVectorizable<T, ?>> extends KernelLib {

	private static final int defaultLutCost = 200;

	private final DFEVectorType<T> m_inputType;
	private final int m_numDimensions;
	private final BoxBufferParams m_1dParams;
	private final List<Buffer1D> buffers;
	private final int[] m_maxItems;
	private final int[] m_numOutputItems;
	private final int[] m_skip;

	private boolean m_hasWritten = false;

	public BoxBuffer(KernelLib root, int maxItems, int numOutputItems, DFEVectorType<T> inputType) {
		this(root, new int[] {maxItems}, new int[] {numOutputItems}, inputType);
	}

	public BoxBuffer(KernelLib root, int maxItems, int numOutputItems, DFEVectorType<T> inputType, boolean doubleBuffered) {
		this(root, new int[] {maxItems}, new int[] {numOutputItems}, inputType, doubleBuffered);
	}

	public BoxBuffer(KernelLib root, int maxItems, int numOutputItems, DFEVectorType<T> inputType, int costOfBramInLuts) {
		this(root, new int[] {maxItems}, new int[] {numOutputItems}, inputType, costOfBramInLuts);
	}

	public BoxBuffer(KernelLib root, int maxItems, int numOutputItems, DFEVectorType<T> inputType, boolean doubleBuffered, int costOfBramInLuts) {
		this(root, new int[] {maxItems}, new int[] {numOutputItems}, inputType, doubleBuffered, costOfBramInLuts);
	}

	public BoxBuffer(KernelLib root, int[] maxItems, int[] numOutputItems, DFEVectorType<T> inputType) {
		this(root, maxItems, numOutputItems, inputType, defaultLutCost);
	}

	public BoxBuffer(KernelLib root, int[] maxItems, int[] numOutputItems, DFEVectorType<T> inputType, boolean doubleBuffered) {
		this(root, maxItems, numOutputItems, inputType, doubleBuffered, defaultLutCost);
	}

	public BoxBuffer(KernelLib root, int[] maxItems, int[] numOutputItems, DFEVectorType<T> inputType, int costOfBramInLuts) {
		this(root, maxItems, numOutputItems, inputType, true, costOfBramInLuts);
	}

	public BoxBuffer(KernelLib root, int[] maxItems, int[] numOutputItems, DFEVectorType<T> inputType, boolean doubleBuffered, int costOfBramInLuts) {
		super(root);
		if (maxItems.length != numOutputItems.length) {
			throw new RuntimeException("The number of dimensions of the box must be the same as the number of dimensions of the buffer.");
		}
		if (maxItems.length > 1) {
			throw new RuntimeException("Only 1D supported right now.");
		}

		m_numDimensions  = maxItems.length;
		m_inputType      = inputType;
		m_numOutputItems = new int[m_numDimensions];
		m_maxItems       = new int[m_numDimensions];
		m_skip           = new int[m_numDimensions];
		for (int i = 0; i < m_numDimensions; i++) {
			m_numOutputItems[i] = numOutputItems[i];
			m_maxItems[i] = maxItems[i];
		}

		int numBuffers = 1;
		int bufferDepth = maxItems[m_numDimensions - 1];
		m_skip[m_numDimensions - 1] = 1;
		for (int i = m_numDimensions - 2; i >= 0; i--) {
			numBuffers  *= numOutputItems[i];
			bufferDepth *= MathUtils.ceilDivide(maxItems[i], numOutputItems[i]);
			m_skip[i]    = MathUtils.ceilDivide(maxItems[i], numOutputItems[i]) * m_skip[i + 1];
		}

		m_1dParams = new BoxBufferParams(bufferDepth, inputType.getSize(), numOutputItems[m_numDimensions - 1], inputType.getContainedType().getTotalBits(), costOfBramInLuts, doubleBuffered);

		buffers = new ArrayList<Buffer1D>(numBuffers);
		for (int i = 0; i < numBuffers; i++) {
			buffers.add(new Buffer1D(inputType.getContainedType()));
		}
	}


	@Override
	protected void finalize() throws Throwable {
	    super.finalize();
	    if (!m_hasWritten) {
	    	throw new RuntimeException("You must write data into the Box Buffer using the write method to use it.");
	    }
	}


	public void write(DFEVector<T> data, DFEVar address, DFEVar enable) {
		write(data, new DFEVar[] {address}, enable);
	}

	public void write(DFEVector<T> data, DFEVar address, DFEVar enable, DFEVar buffer) {
		write(data, new DFEVar[] {address}, enable, buffer);
	}

	public void write(DFEVector<T> data, DFEVar[] address, DFEVar enable) {
		if (m_1dParams.doubleBuffered) {
			throw new RuntimeException("If the Box Buffer is double buffered, then you must specify which buffer you want to write to.");
		}
	}

	public void write(DFEVector<T> data, DFEVar[] address, DFEVar enable, DFEVar buffer) {
		if (!data.getType().equals(m_inputType)) {
			throw new RuntimeException("Type given to write function does not match type passed to the BoxBuffer constructor.");
		}
		if (address.length != m_numDimensions) {
			throw new RuntimeException("You must supply an address for each dimension, i.e. the length of address must be that same as the number of dimensions specified in the constructor.");
		}
		if (m_hasWritten) {
			throw new RuntimeException("The box buffer can only be written to once.");
		}
		m_hasWritten = true;
		optimization.pushEnableBitGrowth(false);//Make sure we don't have any weird modes on

		DFEVar wrRowCol = xyLookup(get1dAddress(address), dfeUInt(m_1dParams.rowAddrBits + m_1dParams.colAddrBits));

		DFEVar wrCol  = slice(wrRowCol, 0, m_1dParams.colAddrBits);
		DFEVar _wrRow = slice(wrRowCol, m_1dParams.colAddrBits, m_1dParams.rowAddrBits);

		// Widen wrRow in if we've had to pad the ram depth
		DFEVar ramWriteCol = wrCol;
		DFEVar ramWriteRow = _wrRow.cast(dfeUInt(m_1dParams.rowAddrBits));
		DFEVar ramWriteRowp1 = ramWriteRow + 1;

		//Stop writing at the end of the RAM
		DFEVar wrEnablep1 = ramWriteRow < m_1dParams.maxAddress ? enable : 0;

		DFEVector<T> paddedInputData = padInput(data);

		DFEVector<T> alignedInputData = alignInput(paddedInputData, ramWriteCol);

		DFEVar[] writeAddress = new DFEVar[m_1dParams.numTiles];
		DFEVar[] writeEnable  = new DFEVar[m_1dParams.numTiles];
		for (int i = 0; i < m_1dParams.numTiles; i++) {
			int tileLastItem = ((i + 1) * m_1dParams.tileItems) - 1;

			// Write port - wrBuffer:ramWriteRow
			DFEVar nextRow = ramWriteCol > tileLastItem;
			DFEVar writeRow = ~nextRow ? ramWriteRow : ramWriteRowp1;
			writeEnable[i] = ~nextRow ? enable : wrEnablep1;

			// Create actual write address
			writeAddress[i] = m_1dParams.doubleBuffered
			                ? writeRow.cat(buffer).cast(dfeUInt(m_1dParams.ramAddrBits))
			                : writeRow.cast(dfeUInt(m_1dParams.ramAddrBits));
		}

		for (int i = 0; i < buffers.size(); i++) {
			DFEVar inThisBuffer = constant.var(true);//TODO: multidim
			buffers[i].write(alignedInputData, writeAddress, writeEnable, inThisBuffer);

		}
		optimization.popEnableBitGrowth();
	}


	public DFEVector<T> read(DFEVar address) {
		return read(new DFEVar[] {address});
	}

	public DFEVector<T> read(DFEVar address, DFEVar buffer) {
		return read(new DFEVar[] {address}, buffer);
	}

	public DFEVector<T> read(DFEVar address[]) {
		if (m_1dParams.doubleBuffered) {
			throw new RuntimeException("If the Box Buffer is double buffered, then you must specify which buffer you want to read from.");
		}
		return read(address, null);
	}

	public DFEVector<T> read(DFEVar[] address, DFEVar buffer) {
		return buffers[0].read(address[0], buffer);//TODO: multidim
	}


	//Slice bits out of a DFEVar and reinterpret as an unsigned integer (rather than raw bits).
	private DFEVar slice(DFEVar input, int base, int width) {
		DFEVar output = input.slice(base, width);
		return output.cast(dfeUInt(width));
	}


	private DFEVar get1dAddress(DFEVar[] address) {
		optimization.pushEnableBitGrowth(true);
		List<DFEVar> summands = new ArrayList<DFEVar>();
		summands.add(address[m_numDimensions - 1]);
		for (int i = 0; i < m_numDimensions - 1; i++) {
			summands.add(address[i] * m_skip[i]);
		}
		DFEVar output = adderTree(summands);
		optimization.popEnableBitGrowth();
		return output.cast(dfeUInt(MathUtils.bitsToAddress(m_1dParams.maxItems)));
	}


	private DFEVar adderTree(List<DFEVar> input) {
		if (input.size() == 1) {
			return input[0];
		}
		List<DFEVar> output = new ArrayList<DFEVar>();
		for (int i = 0; i < input.size() / 2; i++) {
			output.add(input[2 * i] + input[2 * i + 1]);
		}
		if (input.size() % 2 != 0) {
			output.add(input[input.size() - 1]);
		}
		return adderTree(output);
	}


	private DFEVector<T> padInput(DFEVector<T> input) {
		List<T> paddedInputData = new ArrayList<T>();
		for (int i = 0; i < m_1dParams.numCols; i++) {
			if(i < m_1dParams.numInputItems) {
				// New Data
				paddedInputData.add(input[i]);
			} else if(i > (m_1dParams.numCols - m_1dParams.tileItems)) {
				//Old data
				int offsetFromEnd = (m_1dParams.numCols - 1) - i;
				int index = m_1dParams.numInputItems - (1 + (offsetFromEnd % m_1dParams.numInputItems));
				int lookback = 1 + offsetFromEnd / m_1dParams.numInputItems;

				paddedInputData.add(stream.offset(input[index], -lookback));
			} else {
				// Pad
				paddedInputData.add(constant.zero(input[0].getType()));
			}
		}
		return asDFEVector(paddedInputData);
	}


	private DFEVector<T> alignInput(DFEVector<T> input, DFEVar ramWriteCol) {
		DFEVector<T> output;
		if(m_1dParams.numInputItems != m_1dParams.numCols) {
			int gcd = MathUtils.greatestCommonDivisor(m_1dParams.numCols, m_1dParams.numInputItems);
			int chunks = m_1dParams.numCols / gcd;

			DFEVector<DFEVector<T>> chunkedInputData = as2dDFEVector(input, gcd, chunks);
			DFEVar rotate = ConstDenominator.divMod(ramWriteCol, gcd, MathUtils.bitsToAddress(chunks)).m_quot;
			output = asSingleDFEVector(chunkedInputData.rotateElementsLeft(rotate));
		} else {
			output = input;
		}
		return output;
	}


	private DFEVector<T> asDFEVector(List<T> input) {
		DFEVectorType<T> type = new DFEVectorType<T>(input[0].getType(), input.size());
		DFEVector<T> output = type.newInstance(this);
		for (int i = 0; i < input.size(); i++) {
			output[i] <== input[i];
		}
		return output;
	}


	private DFEVector<T> asSingleDFEVector(DFEVector<DFEVector<T>> input) {
		DFEVectorType<T> type = new DFEVectorType<T>(input[0][0].getType(), input.getSize() * input[0].getSize());
		DFEVector<T> output = type.newInstance(this);
		for (int i = 0; i < input.getSize(); i++) {
			for (int j = 0; j < input[i].getSize(); j++) {
				output[i * input[i].getSize() + j] <== input[i][j];
			}
		}
		return output;
	}


	private DFEVector<DFEVector<T>> as2dDFEVector(DFEVector<T> input, int numItems, int numVectors) {
		DFEVectorType<DFEVector<T>> type = new DFEVectorType<DFEVector<T>>(new DFEVectorType<T>(input[0].getType(), numItems), numVectors);
		DFEVector<DFEVector<T>> output = type.newInstance(this);

		for (int i = 0; i < numVectors; i++) {
			for (int j = 0; j < numItems; j++) {
				output[i][j] <== input[i * numItems + j];
			}
		}
		return output;
	}


	/**
	 * Find the row and column for item i in ordered table (row-major)
	 */
	private DFEVar xyLookup(DFEVar index, DFEType addrType) {
		/* If numcols is a power of two we can slice */
		ConstDivModResult dmr = ConstDenominator.divMod(index, m_1dParams.numCols, addrType.getTotalBits() - m_1dParams.colAddrBits);
		return dmr.m_quot.cat(dmr.m_rem.cast(dfeUInt(m_1dParams.colAddrBits))).cast(addrType);
	}

	private class Buffer1D extends KernelLib {
		private final List<Memory<DFEVector<T>>> m_rams;

		public Buffer1D(KernelType<T> itemType) {
			super(BoxBuffer.this);
			DFEVectorType<T> tileType = new DFEVectorType<T>(itemType, m_1dParams.tileItems);
			m_rams = new ArrayList<Memory<DFEVector<T>>>(m_1dParams.numTiles);
			for (int i = 0; i < m_1dParams.numTiles; i++) {
				m_rams.add(mem.alloc(tileType, m_1dParams.ramDepth));
			}
		}


		public void write(DFEVector<T> alignedInputData, DFEVar[] writeAddress, DFEVar[] writeEn, DFEVar enable) {
			for (int i = 0; i < m_1dParams.numTiles; i++) {
				DFEVector<T> tileInputData = slice(alignedInputData, i * m_1dParams.tileItems, m_1dParams.tileItems);
				m_rams[i].write(optimization.limitFanout(writeAddress[i], 2), tileInputData, optimization.limitFanout(writeEn[i] & enable, 2));
			}
		}


		public DFEVector<T> read(DFEVar address, DFEVar buffer) {
			optimization.pushEnableBitGrowth(false);//Make sure we don't have any weird modes on
			DFEVar addr = xyLookup(address, dfeUInt(m_1dParams.readAddrBits));

			// Slice up into mux select + row select
			DFEVar rdShiftAmt = optimization.limitFanout(slice(addr, 0, m_1dParams.colAddrBits), 30);
			DFEVar ramReadRow = slice(addr, m_1dParams.colAddrBits, m_1dParams.rowAddrBits);


			// Read Logic - calculate the next address when we read off the end of the buffer
			DFEVar ramReadRowp1 = ramReadRow + 1;

			List<T> ramOutputData = new ArrayList<T>();
			for (int i = 0; i < m_1dParams.numTiles; i++) {
				// Choose which row to read from this RAM depending on the column
				// the first item is in (use -1 and gte to avoid constant overflow)
				DFEVar readRow = (i + 1) * m_1dParams.tileItems - 1 >= rdShiftAmt ? ramReadRow : ramReadRowp1;

				// Read port - rdBuffer:readRow
				DFEVar readAddress = m_1dParams.doubleBuffered
				                   ? readRow.cat(buffer).cast(dfeUInt(m_1dParams.ramAddrBits))
				                   : readRow.cast(dfeUInt(m_1dParams.ramAddrBits));

				// Read from the RAM and add it to the output
				ramOutputData.addAll(m_rams[i].read(optimization.limitFanout(readAddress, 2)).getElementsAsList());
			}

			optimization.popEnableBitGrowth();
			return alignRamOutput(rdShiftAmt, ramOutputData);
		}


		//Slice bits out of a DFEVar and reinterpret as an unsigned integer (rather than raw bits).
		private DFEVar slice(DFEVar input, int base, int width) {
			DFEVar output = input.slice(base, width);
			return output.cast(dfeUInt(width));
		}

		private DFEVector<T> slice(DFEVector<T> input, int base, int width) {
			List<T> output = new ArrayList<T>();
			for (int j = 0; j < width; j++) {
				output.add(input[base + j]);
			}
			return asDFEVector(output);
		}


		private DFEVector<T> alignRamOutput(DFEVar rdShiftAmt, List<T> ramOut) {
			DFEVector<T> shiftedRamOut = asDFEVector(ramOut).rotateElementsRight(rdShiftAmt);
			DFEVector<T> outItem = slice(shiftedRamOut, 0, m_1dParams.numOutputItems);

			// Minimize unconnected warnings; these will optimize away
			for (int j = outItem.getSize(); j < shiftedRamOut.getSize(); j++) {
				shiftedRamOut[j].setReportOnUnused(false);
			}

			return outItem;
		}
	}


	private class BoxBufferParams {
		final int numInputItems;
		final int numOutputItems;
		final int ramDepth;
		final int maxAddress;
		final int ramAddrBits;
		final int rowAddrBits;
		final int maxItems;

		final boolean doubleBuffered;

		final int colAddrBits;
		final int readAddrBits;
		final int numCols;
		final int numTiles;
		final int tileItems;

		BoxBufferParams(int maxItems, int numInputItems, int numOutputItems, int itemWidthBits, int bramCostInLuts, boolean doubleBuffered) {
			this.doubleBuffered = doubleBuffered;
			this.numInputItems  = numInputItems;
			this.numOutputItems = numOutputItems;
			this.maxItems       = maxItems;

			if(maxItems % numInputItems != 0) {
				throw new MaxCompilerAPIError("maxItems (" + maxItems
					+ ") is not whole multiple of numInputItems ("
					+ numInputItems + ")");
			}

			int numBuffers = doubleBuffered ?  2 :  1;
			tileItems = getTileWidth(itemWidthBits, numBuffers, maxItems, numInputItems, numOutputItems, bramCostInLuts);

			// This is the total width needed
			int totalWidthNeeded = getTotalWidthNeeded(tileItems, numInputItems, numOutputItems);
			numTiles = MathUtils.ceilDivide(totalWidthNeeded, tileItems);

			// Assume all the tiles have the same number of items
			numCols     = numTiles * tileItems;
			colAddrBits = MathUtils.bitsToAddress(numCols);

			// Number of rows, such that rows*cols >= maxItems and so that it lines up with a multiple of numInputItems
			int paddedItems = MathUtils.nextMultiple(maxItems, MathUtils.leastCommonMultiple(numCols, numInputItems));
			maxAddress  = MathUtils.ceilDivide(paddedItems, numCols) - 1;
			ramDepth    = numBuffers * (maxAddress+ 1);
			ramAddrBits = MathUtils.bitsToAddress(ramDepth);
			rowAddrBits = ramAddrBits - MathUtils.bitsToAddress(numBuffers);

			readAddrBits = rowAddrBits + colAddrBits;
		}

		int estimateMinRamWidth(long depth, int bramWidth) {
			int numOptions = MathUtils.bitsToRepresent(MathUtils.greatestCommonDivisor(bramWidth, MathUtils.nextPowerOfTwo(bramWidth)));
			for (int i = 0; i < numOptions; i++) {
				long optionDepth = 512       << (numOptions - 1 - i);
				int  optionWidth = bramWidth >> (numOptions - 1 - i);
				if (depth % optionDepth == 0) {
					return optionWidth;
				}
			}
			return bramWidth;
		}

		int getTotalWidthNeeded(int tileItems, int numInputItems, int numOutputItems) {
			// What do the input / outputs need?
			int inputNeeded  = numInputItems  + Math.min(numInputItems  - 1, tileItems - 1);
			int outputNeeded = numOutputItems + Math.min(numOutputItems - 1, tileItems - 1);

			// This is the total width needed
			return Math.max(inputNeeded, outputNeeded);
		}

		int getTileWidth(int itemWidthBits, int numBuffers, int maxItems, int numInputItems, int numOutputItems, int costOfBramInLuts) {
			int bestWidth = -1;
			int bestScore= Integer.MAX_VALUE;

			boolean isAltera = BoxBuffer.this.getManager().getManagerConfiguration().getBoardModel().isAltera();
			int bramWidth = isAltera ? 40 : 72;

			int biggestTile = MathUtils.leastCommonMultiple(bramWidth, itemWidthBits);
			int mostItems = biggestTile / itemWidthBits;

			for (int i = 1; i <= mostItems; i++) {
				int bramCost = costInBrams(i, itemWidthBits, numBuffers, maxItems, numInputItems, numOutputItems, bramWidth);
				int lutCost  = costInLuts(i, itemWidthBits, numInputItems, numOutputItems, maxItems);
				int score = lutCost + costOfBramInLuts * bramCost;
				if(score < bestScore) {
					bestScore = score;
					bestWidth = i;
				}
			}
			return bestWidth;
		}

		int costInBrams(int tileItems, int itemWidthBits, int numBuffers, int maxItems, int numInputItems, int numOutputItems, int bramWidth) {
			int totalWidthNeeded = getTotalWidthNeeded(tileItems, numInputItems, numOutputItems);
			int depth = MathUtils.nextMultiple(numBuffers * MathUtils.ceilDivide(maxItems, totalWidthNeeded), 512);
			int optWidth = estimateMinRamWidth(depth, bramWidth);

			int tileWidth = MathUtils.nextMultiple(itemWidthBits * tileItems, optWidth);
			int numTiles = MathUtils.ceilDivide(totalWidthNeeded, tileItems);

			return MathUtils.ceilDivide(depth * tileWidth * numTiles, bramWidth * 512);
		}

		//TODO: move everything to the finaliser so that we know how many times it has been read from and multiply the number of LUTs for output by that.
		int costInLuts(int tileItems, int itemWidthBits, int numInputItems, int numOutputItems, int maxItems) {
			int totalWidthNeeded = getTotalWidthNeeded(tileItems, numInputItems, numOutputItems);
			int numTiles = MathUtils.ceilDivide(totalWidthNeeded, tileItems);

			int numCols = numTiles * tileItems;

			int gcd = MathUtils.greatestCommonDivisor(numCols, numInputItems);
			int chunks = numCols / gcd;
			int bitsForNumerator = MathUtils.bitsToAddress(maxItems / numInputItems);
			int bitsForDenominator = MathUtils.bitsToRepresent(gcd);

			int lutsForOutput = getLutsForRotate(numCols, itemWidthBits);
			int lutsForInput  = getLutsForRotate(chunks, itemWidthBits * gcd);
			lutsForInput += MathUtils.isPowerOf2(gcd) ? 0 : bitsForNumerator * bitsForDenominator;//LUTs for DivMod
			lutsForInput = numCols == numInputItems ? 0 : lutsForInput;//No rotate needed in this case.

			return lutsForOutput + lutsForInput;
		}

		int getLutsForRotate(int size, int bitWidth) {
			boolean isAltera = BoxBuffer.this.getManager().getManagerConfiguration().getBoardModel().isAltera();
			int luts;
			if (isAltera) {
				luts = Math.min(log(size, 4) * 2, log(size, 3));
			} else {
				luts = log(size, 4);
			}

			return luts * size * bitWidth;
		}

		int log(int x, int base) {
			return (int)Math.ceil(Math.log(x) / Math.log(base));
		}
	}
}
