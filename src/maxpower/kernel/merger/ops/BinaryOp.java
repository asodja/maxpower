package maxpower.kernel.merger.ops;

/**
 * Binary operation.
 *
 * @author nvoss
 */
public interface BinaryOp {
	/**
	 * Evaluates the operation on two constants.
	 * @param val1 First input value.
	 * @param val2 Second input value.
	 * @return Result.
	 */
	public double evaluate(double val1, double val2);
}
