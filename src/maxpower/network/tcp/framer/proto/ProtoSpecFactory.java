/*********************************************************************
 * TCP Framer                                                        *
 * Copyright (C) 2013-2015 Maxeler Technologies                      *
 *                                                                   *
 * Author:  Itay Greenspon                                           *
 *                                                                   *
 *********************************************************************/

package maxpower.network.tcp.framer.proto;

import maxpower.network.tcp.framer.TCPFramerSM;

public interface ProtoSpecFactory {
	public FramerProtocolSpec create(TCPFramerSM owner);
}
