import { Line } from '@react-three/drei'
import { useMemo } from 'react'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO, LOOP_RADIUS_KM, EARTH_RADIUS_KM } from '../../config/simulation'

const LoopTethers = () => {
  const tetherLines = useMemo(() => {
    return new Array(6).fill(0).map((_, idx) => {
      const t = idx / 5
      const anchor = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t, LOOP_RADIUS_KM)
      const ground = anchor.clone().setLength(EARTH_RADIUS_KM)
      return [anchor, ground]
    })
  }, [])

  return (
    <group>
      {tetherLines.map((pair, idx) => (
        <Line key={idx} points={pair.map((p) => p.toArray())} color="#b8c0ff" lineWidth={1.5} />
      ))}
    </group>
  )
}

export default LoopTethers
