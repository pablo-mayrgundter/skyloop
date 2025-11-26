import {Line} from '@react-three/drei'
import {useMemo} from 'react'
import {Vector3} from 'three'
import {LOS_ANGELES, SAN_FRANCISCO} from '../../config/simulation'
import {buildGreatCircleArc} from '../Globe'

type LoopTethersProps = {
  count?: number
  lengthKm?: number
}

const LoopTethers = ({ count = 12, lengthKm = 20 }: LoopTethersProps) => {
  const tetherSegments = useMemo(() => {
    const arc = buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, count)
    return arc.map((p) => ({
      start: p.clone().add(new Vector3(0, lengthKm * 0.001, 0)),
      end: p
    }))
  }, [count, lengthKm])

  return (
    <group>
      {tetherSegments.map((segment, idx) => (
        <Line
          key={idx}
          points={[segment.start, segment.end]}
          color="#9be0a9"
          lineWidth={1}
          transparent
          opacity={0.6}
        />
      ))}
    </group>
  )
}

export default LoopTethers
