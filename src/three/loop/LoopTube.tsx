import { Line } from '@react-three/drei'
import { useMemo } from 'react'
import { LOS_ANGELES, SAN_FRANCISCO } from '../../config/simulation'
import { buildGreatCircleArc } from '../Globe'

type LoopTubeProps = {
  segments?: number
}

const LoopTube = ({ segments = 64 }: LoopTubeProps) => {
  const tubePoints = useMemo(() => buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, segments), [segments])
  return <Line points={tubePoints} color="#ffa3c4" lineWidth={3} />
}

export default LoopTube
