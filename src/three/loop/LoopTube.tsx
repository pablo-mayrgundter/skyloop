import { Tube } from '@react-three/drei'
import { useMemo } from 'react'
import { CatmullRomCurve3, Vector3 } from 'three'
import { LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle, LOOP_RADIUS_KM } from '../../config/simulation'

const LoopTube = () => {
  const curve = useMemo(() => {
    const points: Vector3[] = []
    const segments = 24
    for (let i = 0; i <= segments; i += 1) {
      points.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, i / segments, LOOP_RADIUS_KM))
    }
    return new CatmullRomCurve3(points)
  }, [])

  // Tube radius: 0.0015 km = 1.5 meters (reasonable for a vacuum tube)
  return (
    <Tube args={[curve, 128, 0.0015, 12, false]}>
      <meshStandardMaterial color="#6aa0ff" transparent opacity={0.4} />
    </Tube>
  )
}

export default LoopTube
