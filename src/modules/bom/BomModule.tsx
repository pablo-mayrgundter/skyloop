import { CostBreakdown } from '../../types/simulation'

const items: CostBreakdown[] = [
  { category: 'Sections', costMillions: 4200 },
  { category: 'Tethers', costMillions: 1100 },
  { category: 'Power', costMillions: 800 },
  { category: 'Stations', costMillions: 650 }
]

const BomModule = () => {
  const total = items.reduce((sum, item) => sum + item.costMillions, 0)
  return (
    <div className="controls-grid">
      <h3>Bill of Materials</h3>
      {items.map((item) => (
        <div key={item.category} className="metric">
          <strong>{item.category}</strong>
          <div>${item.costMillions.toLocaleString()}M</div>
        </div>
      ))}
      <div className="metric">
        <strong>Total</strong>
        <div>${total.toLocaleString()}M</div>
      </div>
      <p>TODO: connect to parametric cost estimators and sensitivity analyses.</p>
    </div>
  )
}

export default BomModule
