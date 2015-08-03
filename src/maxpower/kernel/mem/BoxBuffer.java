package maxpower.kernel.mem;

import com.maxeler.maxcompiler.v2.errors.MaxCompilerAPIError;
import com.maxeler.maxcompiler.v2.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v2.kernelcompiler.KernelLib;
import com.maxeler.maxcompiler.v2.kernelcompiler._KernelConfiguration;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.KernelMath;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.KernelMath.DivModResult;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem.RamPortMode;
import com.maxeler.maxcompiler.v2.kernelcompiler.stdlib.core.Mem.RamWriteMode;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.KernelObjectVectorizable;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.KernelType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEType;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.base.DFEVar;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVector;
import com.maxeler.maxcompiler.v2.kernelcompiler.types.composite.DFEVectorType;
import com.maxeler.maxcompiler.v2.utils.MathUtils;
import com.maxeler.photon.configuration.PhotonPlatform;

/**
 * The Box buffer reads in a certain sized stream of data and holds it in a
 * buffer, from which you can query it to retrieve a box of adjacent items at
 * a random location. In this case box is meant in the geometric sense (also
 * called n-orthotope or hyperectangle), so a contiguous section of an N-dimensional
 * space can be retrieved. In 1D this is a simple interval.
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
 */
public class BoxBuffer<T extends KernelObjectVectorizable<T, ?>> extends KernelLib {

	private final Kernel root;
	static int xyIndex = 0;

	public static class BoxBufferParams {
		int maxItems;
		final int numBuffers;
		final static int minRamDepth = 512;
		final public int itemWidthBits;
		final public int numInputItems;
		final public int numOutputItems;
		final int ramDepth;
		final public int ramAddrBits;
		final public int rowAddrBits;
		final public int bufAddrBits;

		final public int colAddrBits;

		final public int numRows;
		final public int numCols;

		final public int numTiles;
		final public int tileItems;

		final private Kernel root;

		/** Number of bits for the read address to the WDR */
		final public int readAddrBits;

		public BoxBufferParams(int maxItems, int numInputItems,
			int numOutputItems, int itemWidthBits, Kernel root) {
			this(maxItems, numInputItems, numOutputItems, itemWidthBits, 2, root);
		}

		public BoxBufferParams(int maxItems, int numInputItems,
			int numOutputItems, int itemWidthBits, int numBuffers, Kernel root) {

			this.numBuffers = numBuffers;
			this.itemWidthBits = itemWidthBits;
			this.numInputItems = numInputItems;
			this.numOutputItems = numOutputItems;
			this.maxItems = maxItems;
			this.root = root;

			if(maxItems % numInputItems != 0)
				throw new MaxCompilerAPIError("maxItems (" + maxItems
					+ ") is not whole multiple of numInputItems ("
					+ numInputItems + ")");
			tileItems = getTileWidth();

			// This is the total width needed
			int totalWidthNeeded = getTotalWidthNeeded(tileItems);

			numTiles = (int) Math.ceil((double) totalWidthNeeded / tileItems);

			// Assume all the tiles have the same number of items
			numCols = numTiles * tileItems;
			colAddrBits = MathUtils.bitsToAddress(numCols);

			// Number of rows, such that rows*cols >= maxItems
			numRows = MathUtils
				.nextPowerOfTwo((int) ((double) maxItems / numCols));
			ramDepth = Math.max(minRamDepth, MathUtils
				.nextPowerOfTwo(numBuffers * maxItems / numCols));
			ramAddrBits = MathUtils.bitsToAddress(ramDepth);
			bufAddrBits = MathUtils.bitsToAddress(numBuffers);
			rowAddrBits = ramAddrBits - bufAddrBits;

			readAddrBits = rowAddrBits + colAddrBits;
		}

		int estimateDepth(int minDepth) {
			return (int) Math.ceil(Math.pow(2.0, Math.max(9.0, Math.ceil(Math
				.log(minDepth)
				/ Math.log(2.0)))));
		}

		int estimateMinRamWidth(int depth) {
			int bramWidth;
			if (_KernelConfiguration.getPhotonKernelConfig(root.getManager().getCurrentKernelConfig()).getPlatform()
				== PhotonPlatform.AlteraDefault)
				bramWidth = 40;
			else
				bramWidth = 36; //TODO: I think this should be 72
			return (int) Math.floor(bramWidth / Math.pow(2, Math.max(0, (Math
				.log(depth) / Math.log(2)) - 9.0)));
		}

		int getTotalWidthNeeded(int tileItems) {
			// What do the input / outputs need?
			int inputNeeded = numInputItems
				+ Math.min(this.numInputItems - 1, tileItems - 1);
			int outputNeeded = numOutputItems
				+ Math.min(this.numOutputItems - 1, tileItems - 1);

			// This is the total width needed
			return Math.max(inputNeeded, outputNeeded);
		}

		int getTileWidth() {
			int bestWidth = -1;
			double bestScore = 0.0;
			int optimalBramWidth;
			if (_KernelConfiguration.getPhotonKernelConfig(root.getManager().getCurrentKernelConfig()).getPlatform()
				== PhotonPlatform.AlteraDefault)
				optimalBramWidth = 40;
			else
				optimalBramWidth = 72;

			int biggestTile = MathUtils.leastCommonMultiple(optimalBramWidth, itemWidthBits);
			int mostItems = biggestTile / itemWidthBits;

			for (int i = 1; i <= mostItems; i++) {
				double score = estimateEfficiency(i);
				if(score > bestScore) {
					bestScore = score;
					bestWidth = i;
				}
			}
			return bestWidth;
		}

		double estimateEfficiency(int tileItems) {
			int totalWidthNeeded = getTotalWidthNeeded(tileItems);
			int depth = estimateDepth((numBuffers * maxItems + (totalWidthNeeded / 2))
				/ totalWidthNeeded);
			int optWidth = estimateMinRamWidth(depth);

			double packingEff = (double) itemWidthBits
				* Math.max(numInputItems, numOutputItems)
				/ (MathUtils.nextMultiple(itemWidthBits * tileItems, optWidth)
					* MathUtils.ceilDivide(totalWidthNeeded, tileItems));

			return packingEff;
		}

	}

	public BoxBuffer(Kernel _root) {
		super(_root);
		root = _root;
	}

	public DFEVector<T> makeTiledRam(
		BoxBufferParams params,
		DFEVector<T> wrData,
		DFEVar wrBuffer,
		DFEVar wrRow,
		DFEVar wrCol,
		DFEVar wrEnable,
		DFEVar rdBuffer,
		DFEVar ramReadRow,
		DFEVar rdColSelect)
	{
		/*
		 * Striped version of the RAM
		 *
		 * Each input item is stored in a separate ram, we write [numInputItems]
		 * at a time, we expect to have a whole-multiple of numInputItems RAMs
		 */

		KernelType<T> _itemType = wrData[0].getType();
		DFEVectorType<T> arrayType = new DFEVectorType<T>(_itemType, params.numCols);
		DFEVector<T> ramOutputData = arrayType.newInstance(root);

		// Widen wrRow in if we've had to pad the ram depth
		DFEVar ramWriteCol = wrCol;
		DFEVar ramWriteRow = wrRow.cast(dfeUInt(params.rowAddrBits));// ;
		DFEVar ramWriteRowp1 = ramWriteRow + 1;

		// Make wrEnable, prevent next rows writing at 0 if we're right at the
		// end of the RAM
		DFEVar wrEnablep1 = ramWriteRowp1.gt(constant.var(0)) ? wrEnable
			: constant.zero(dfeBool());

		// Read Logic - calculate the next address when we read off the end of
		// the buffer
		DFEVar ramReadRowp1 = ramReadRow + 1;

		DFEVector<T> paddedInputData = arrayType.newInstance(root);

		/**
		 * Build up the input data word we use to write into the RAMs
		 *
		 * ----------------------------------------- | 0 | 1 | 2 | 3 | X | X | X
		 * |-3 |-2 |-1 | ----------------------------------------- | new data |
		 * pad | old data | -----------------------------------------
		 *
		 * We put the new data in, plus tileWidth-1 items of oldData
		 *
		 */

		for (int i = 0; i < params.numCols; i++) {
			if(i < params.numInputItems) {
				// New Data
				paddedInputData[i] <== wrData[i];
			} else if(i > (params.numCols - params.tileItems)) {

				int offsetFromEnd = (params.numCols - 1) - i;
				int index = params.numInputItems
					- (1 + (offsetFromEnd % params.numInputItems));
				int lookback = 1 + offsetFromEnd / params.numInputItems;

				paddedInputData[i] <== stream.offset(wrData[index], -lookback);
			} else {
				// Pad
				paddedInputData[i] <== constant.zero(_itemType);
			}
		}

		DFEVector<T> alignedInputData = null;

		if(params.numInputItems != params.numCols) {
			int gcd = MathUtils.greatestCommonDivisor(params.numCols,
				params.numInputItems);
			int chunks = params.numCols / gcd;

			DFEVectorType<DFEVector<T>> chunkedType = new DFEVectorType<DFEVector<T>>(
				new DFEVectorType<T>(_itemType, gcd), chunks);
			DFEVector<DFEVector<T>> chunkedInputData = chunkedType.newInstance(root);

			for (int i = 0; i < chunks; i++)
				for (int j = 0; j < gcd; j++)
					chunkedInputData[i][j] <== paddedInputData[i * gcd + j];

			DivModResult dmr = KernelMath.divMod(ramWriteCol,
				constant.var(gcd), MathUtils.bitsToAddress(chunks));
			DFEVar rotate = dmr.getQuotient();

			DFEVector<DFEVector<T>> rotatedInputData = chunkedInputData
				.rotateElementsLeft(rotate);
			alignedInputData = arrayType.newInstance(root);

			for (int i = 0; i < chunks; i++)
				for (int j = 0; j < gcd; j++)
					alignedInputData[i * gcd + j] <== rotatedInputData[i][j];
		} else
			alignedInputData = paddedInputData;

		for (int i = 0; i < params.numTiles; i++) {

			int tileLastItem = ((i + 1) * params.tileItems) - 1;

			// Write port - wrBuffer:ramWriteRow
			DFEVar nextRow = ramWriteCol > tileLastItem;
			DFEVar writeRow = ~nextRow ? ramWriteRow : ramWriteRowp1;
			DFEVar writeEn = ~nextRow ? wrEnable : wrEnablep1;

			DFEVectorType<T> tileType = new DFEVectorType<T>(_itemType,
				params.tileItems);
			DFEVector<T> tileInputData = tileType.newInstance(root);

			// Connect up data for this tile
			int tileItemCol = i * params.tileItems;
			for (int j = 0; j < params.tileItems; j++)
				tileInputData[j] <== alignedInputData[tileItemCol + j];

			// Create actual write address
			DFEVar writeAddress =
				(params.numBuffers > 1) ?
				wrBuffer.cat(writeRow).cast(dfeUInt(params.ramAddrBits))
				: writeRow.cast(dfeUInt(params.ramAddrBits));

			Mem.RamPortParams<DFEVector<T>> portA = mem.makeRamPortParams(
				RamPortMode.WRITE_ONLY,
				optimization.limitFanout(writeAddress, 2), tileType)
				.withDataIn(tileInputData).withWriteEnable(
					optimization.limitFanout(writeEn, 2));

			// Choose which row to read from this RAM depending on the column
			// the first item is in (use -1 and gte to avoid constant overflow)
			DFEVar readRow = (i + 1) * params.tileItems - 1 >= rdColSelect ? ramReadRow
				: ramReadRowp1;

			// Read port - rdBuffer:readRow
			DFEVar readAddress =
				(params.numBuffers > 1) ?
				rdBuffer.cat(readRow).cast(dfeUInt(params.ramAddrBits)):
					readRow.cast(dfeUInt(params.ramAddrBits));
			Mem.RamPortParams<DFEVector<T>> portB = mem.makeRamPortParams(
				RamPortMode.READ_ONLY,
				optimization.limitFanout(readAddress, 2), tileType);

			// Create the RAM and add it to the output array
			DFEVector<T> tileOutputData = mem.ramDualPort(params.ramDepth,
				RamWriteMode.READ_FIRST, portA, portB).getOutputB();

			for (int j = 0; j < params.tileItems; j++) {
				int outputItem = i * params.tileItems + j;
				ramOutputData[outputItem] <== tileOutputData[j];
			}

		}

		return ramOutputData;

	}

	public DFEVector<T> alignRamOutput(
		BoxBufferParams params,
		DFEVar rdShiftAmt,
		DFEVector<T> ramOut,
		DFEVectorType<T> outType)
	{
		DFEVector<T> shiftedRamOut = ramOut.rotateElementsRight(rdShiftAmt);
		DFEVector<T> outItem = outType.newInstance(root);
		for (int i = 0; i < outItem.getSize(); i++)
			outItem[i] <== shiftedRamOut[i];

		// Minimize unconnected warnings; these will optimize away
		for (int j = outItem.getSize(); j < ramOut.getSize(); j++)
			shiftedRamOut[j].setReportOnUnused(false);

		return outItem;
	}

	/**
	 * Find the row and column for item i in ordered table (row-major)
	 */
	DFEVar xyLookup(
		DFEVar rdIndex,
		DFEType addrType,
		int numRows,
		int numCols,
		int colAddrBits)
	{
		/* If numcols is a power of two we can slice */
		if(MathUtils.isPowerOf2(numCols)) {
			if((int) Math.pow(2, colAddrBits) != numCols)
				throw new MaxCompilerAPIError(getManager(), "Math.pow(2, colAddrBits) != numCols");
			else
				return rdIndex.cast(addrType);
		} else {
			DivModResult dmr = KernelMath.divMod(rdIndex,
				constant.var(numCols), addrType.getTotalBits() - colAddrBits);
			return dmr.getQuotient().cat(
				dmr.getRemainder().cast(dfeUInt(colAddrBits))).cast(addrType);
		}
	}

	/**
	 * We have a stream N items wide, coming into this block
	 *
	 * We want to be able to read a window of M items from this stream out at
	 * any point (provided the data has been written to the buffer)
	 *
	 * Double buffering is achieved simply with wrBuffer + rdBuffer bits
	 * indicating which buffer we're writing and reading to.
	 *
	 * wrRow is the row (N items wide) we're writing to, we assume item's 0 to
	 * N-1 are in row 0, N to 2N-1 in row 1, etc.
	 *
	 * wrData is the N item wide data stream
	 *
	 * rdIndex is the index of the first of M consecutive values we want to read
	 *
	 * wrEnable is the write-enable for writing to the buffer
	 */
	public DFEVector<T> readBox(
		DFEVar wrBuffer,
		DFEVar wrRow,
		DFEVector<T> wrData,
		DFEVar rdBuffer,
		DFEVar rdIndex,
		BoxBufferParams params,
		DFEVar wrEnable)
	{
		// Check bitwidths are okay
		int itemWidthBits = wrData[0].getType().getTotalBits();
		for (int i = 1; i < wrData.getSize(); i++)
			if(wrData[i].getType().getTotalBits() != itemWidthBits)
				throw new MaxCompilerAPIError(
					"All inputs to a window read RAM must be the same width.");

		/*
		 * The address type is..
		 */
		DFEType addrType = dfeUInt(params.readAddrBits);
		DFEVar addr = xyLookup(rdIndex, addrType, params.numRows,
			params.numCols, params.colAddrBits);

		// Slice up into mux select + row select
		DFEVar tmp = addr.slice(0, params.colAddrBits).cast(
			dfeUInt(params.colAddrBits));
		DFEVar rdShiftAmt = optimization.limitFanout(tmp, 30);
		DFEVar ramReadRow = addr.slice(params.colAddrBits, params.rowAddrBits)
			.cast(dfeUInt(params.rowAddrBits));

		KernelType<T> _itemType = wrData[0].getType();
		DFEVectorType<T> outType = new DFEVectorType<T>(_itemType,
			params.numOutputItems);

		DFEType wrRowType = dfeUInt(params.rowAddrBits + params.colAddrBits);
		DFEType wrRow_m_numInputItemsType = dfeUInt(wrRow.getType()
			.getTotalBits()
			+ MathUtils.bitsToRepresent(params.numInputItems));
		DFEVar wrRow_m_numInputItems = wrRow.cast(wrRow_m_numInputItemsType)
			* params.numInputItems;
		DFEVar wrRowCol = xyLookup(wrRow_m_numInputItems, wrRowType,
			params.numRows, params.numCols, params.colAddrBits);

		DFEVar _wrRow = wrRowCol.slice(params.colAddrBits, params.rowAddrBits)
			.cast(dfeUInt(params.rowAddrBits));
		DFEVar wrCol = wrRowCol.slice(0, params.colAddrBits).cast(
			dfeUInt(params.colAddrBits));
		DFEVector<T> ramOut = makeTiledRam(params, wrData, wrBuffer, _wrRow,
			wrCol, wrEnable, rdBuffer, ramReadRow, rdShiftAmt);

		DFEVector<T> outItem = alignRamOutput(params, rdShiftAmt, ramOut, outType);

		return outItem;

	}
}
