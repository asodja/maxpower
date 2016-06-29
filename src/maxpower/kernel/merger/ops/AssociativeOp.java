package maxpower.kernel.merger.ops;

/**
 * Associative binary operation.
 *
 * @author nvoss
 */
public interface AssociativeOp extends BinaryOp {
	/**
	 * Returns the identity of the associative operation. E.g. 0 for addition.
	 * @return The identity of the operation.
	 */
	public double getIdentity();
}
