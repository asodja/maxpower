/*********************************************************************
 * TCP Framer                                                        *
 * Copyright (C) 2013-2015 Maxeler Technologies                      *
 *                                                                   *
 * Author:  Itay Greenspon                                           *
 *                                                                   *
 *********************************************************************/

package maxpower.network.tcp.framer;

import java.util.List;

import maxpower.network.tcp.framer.TcpFramerSM.FramerStates;
import maxpower.network.tcp.framer.proto.FramerProtocolSpec;

import com.maxeler.maxcompiler.v2.managers.DFEManager;
import com.maxeler.maxcompiler.v2.utils.MathUtils;

public class TcpFramerConstants {
	public static boolean enableDebug = false;
	public static boolean enableDebugPrints = false;
	public static boolean enableDebugStreams = false;
	public static boolean enableFramerDebugger = false;

	public static enum FramerErrorCodes {
		NoError,
		HeaderCorrupt,
		PayloadError,

		ShutdownDrain,
		BodyLengthTooBig,

		Reserved3,
		PayloadCutShort,

		PreviousErrors
	}

	public static final int maxSupportedMessageLength = 16 * 1024;
	public static final int modWidth = MathUtils.bitsToAddress(TcpInterfaceTypes.dataWordSizeBytes);
	public static final int outputBufferDepth = MathUtils.nextPowerOfTwo(2 * TcpFramerConstants.maxSupportedMessageLength / TcpInterfaceTypes.dataWordSizeBytes);
	public static final int outputBufferProgrammableFull = outputBufferDepth - (TcpFramerConstants.maxSupportedMessageLength / TcpInterfaceTypes.dataWordSizeBytes + 64 /* fudge factor */);
	public static final int outputEmptyLatency = 16;

	public static void addMaxfileConstants(DFEManager owner, List<FramerProtocolSpec> protoSpecs) {
		owner.addMaxFileConstant("framer_maxSupportedMessageLength", TcpFramerConstants.maxSupportedMessageLength);
		owner.addMaxFileConstant("framer_maxWindowMemorySizeBytes", TcpInterfaceTypes.maxWindowMemorySizeBytes);
		owner.addMaxFileConstant("framer_dataWordSizeBytes", TcpInterfaceTypes.dataWordSizeBytes);
		owner.addMaxFileConstant("framer_outputBufferDepth", TcpFramerConstants.outputBufferDepth);
		owner.addMaxFileConstant("framer_outputBufferProgrammableFull", TcpFramerConstants.outputBufferProgrammableFull);
		owner.addMaxFileConstant("framer_outputEmptyLatency", TcpFramerConstants.outputEmptyLatency);

		owner.addMaxFileConstant("framer_FixFramerConstants_enableDebugStreams", TcpFramerConstants.enableDebugStreams ? 1 : 0);

		owner.addMaxFileConstant("framer_FixFramerConstants_enableDebug", TcpFramerConstants.enableDebug ? 1 : 0);
		owner.addMaxFileConstant("framer_FixFramerConstants_enableFramerDebugger", TcpFramerConstants.enableFramerDebugger ? 1 : 0);

		for (FramerStates state : FramerStates.values()) {
			owner.addMaxFileConstant("framer_FramerStates_" + state.name(), state.ordinal());
		}

		for (FramerErrorCodes errorCode : FramerErrorCodes.values()) {
			owner.addMaxFileConstant("framer_FramerErrorCodes_" + errorCode.name(), errorCode.ordinal());
		}

		for (FramerProtocolSpec pSpec : protoSpecs) {
			owner.addMaxFileConstant("framer_ProtocolSpec_" + pSpec.getProtocolName(), protoSpecs.indexOf(pSpec));
			pSpec.addMaxfileConstants(owner);
		}

	}
}

